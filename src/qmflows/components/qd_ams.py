__all__ = ['check_sys_var', 'ams_job_mopac_sp', 'ams_job_uff_opt']

import os
import shutil

from scm.plams import (Settings, AMSJob, init, finish, Units)
from scm.plams.core.jobrunner import JobRunner
from scm.plams.tools.kftools import KFFile
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)

import scm.plams.interfaces.molecule.rdkit as molkit

from .qd_import_export import export_mol
from .qd_functions import (adf_connectivity, fix_h)


def check_sys_var():
    """
    Check if all ADF environment variables are set and if the 2018 version of ADF is installed.
    """
    sys_var = ['ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            print('WARNING: The environment variable ' + sys_var[i] + ' has not been set')
    if False in sys_var_exists:
        raise EnvironmentError('One or more ADF environment variables have not been set, aborting '
                               'ADF job.')
    if '2018' not in os.environ['ADFHOME']:
        error = 'No ADF version 2018 detected in', os.environ['ADFHOME'], ', aborting ADF job.'
        raise ImportError(error)


def qd_opt(plams_mol, database):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.

    plams_mol <plams.Molecule> The input quantom dot with the 'name' property.
    database <pd.DataFrame>: A database of previous calculations.

    return <plams.Molecule>: An optimized quantom dot.
    """
    name = plams_mol.properties.name.rsplit('.', 1)[0]
    if database is None:
        plams_mol = ams_job_uff_opt(plams_mol)
        plams_mol.properties.entry = False
    elif database.empty or name not in list(database['Quantum_dot_name']):
        plams_mol = ams_job_uff_opt(plams_mol)
        plams_mol.properties.entry = True
    else:
        index = list(database['Quantum_dot_name']).index(name)
        plams_mol_new = molkit.readpdb(database['Quantum_dot_opt_pdb'][index])
        plams_mol_new.properties = plams_mol.properties
        plams_mol = plams_mol_new

    return plams_mol


class CRSResults(SCMResults):
    """
    A class accessing results of COSMO-RS jobs.
    """
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}


class CRSJob(SCMJob):
    """
    A class for running COSMO-RS jobs.
    """
    _command = 'crs'
    _result_type = CRSResults


def uff_settings(plams_mol, mol_indices):
    """
    UFF settings for a constrained geometry optimization
    """
    s1 = Settings()

    s1.input.ams.Task = 'GeometryOptimization'
    s1.input.ams.Constraints.Atom = mol_indices
    s1.input.ams.System.BondOrders._1 = adf_connectivity(plams_mol)
    s1.input.ams.GeometryOptimization.MaxIterations = 100
    s1.input.uff.Library = 'UFF'

    return s1


def crs_settings():
    """
    COSMO-RS settings for a LogP calculation using 'MOPAC PM6' parameters.
    """
    s2 = Settings()
    s2.ignore_molecule = True
    s2.pickle = False

    s2.input.Property._h = 'logP'
    s2.input.Property.VolumeQuotient = 6.766

    s2.input.CRSParameters._1 = 'HB_HNOF'
    s2.input.CRSParameters._2 = 'HB_TEMP'
    s2.input.CRSParameters._3 = 'FAST'
    s2.input.CRSParameters._4 = 'COMBI2005'
    s2.input.CRSParameters.rav = 0.400
    s2.input.CRSParameters.aprime = 1510.0
    s2.input.CRSParameters.fcorr = 2.802
    s2.input.CRSParameters.chb = 8850.0
    s2.input.CRSParameters.sigmahbond = 0.00854
    s2.input.CRSParameters.aeff = 6.94
    s2.input.CRSParameters.Lambda = 0.130
    s2.input.CRSParameters.omega = -0.212
    s2.input.CRSParameters.eta = -9.65
    s2.input.CRSParameters.chortf = 0.816

    s2.input.Compound._h = ''

    s2.input.compound._h = '"' + os.path.join(os.environ['ADFRESOURCES'],
                                              'ADFCRS/1-Octanol.coskf') + '"'
    s2.input.compound.frac1 = 0.725
    s2.input.compound.pvap = 1.01325
    s2.input.compound.tvap = 468.0
    s2.input.compound.meltingpoint = 257.0
    s2.input.compound.flashpoint = 354.0
    s2.input.compound.density = 0.824

    s2.input.compounD._h = '"' + os.path.join(os.environ['ADFRESOURCES'],
                                              'ADFCRS/Water.coskf') + '"'
    s2.input.compounD.frac1 = 0.275
    s2.input.compounD.frac2 = 1.0
    s2.input.compounD.pvap = 1.01325
    s2.input.compounD.tvap = 373.15
    s2.input.compounD.meltingpoint = 273.15
    s2.input.compounD.hfusion = 1.436
    s2.input.compounD.density = 0.997

    return s2


def mopac_settings(charge):
    """
    COSMO-MOPAC settings for a single point with COSMO-RS parameters.
    """
    s3 = Settings()
    s3.pickle = False

    s3.input.ams.Task = 'SinglePoint'
    s3.input.ams.System.Charge = charge
    s3.input.MOPAC.Solvation = 'COSMO-CRS'
    s3.input.MOPAC.mozyme = True

    return s3


def ams_job_mopac_sp(mol_list):
    """
    Runs a MOPAC + COSMO-RS single point.

    mol_list <list>[<plams.Molecule>]: A list of PLAMS molecule with the 'path', 'charge' and
        'name' properties.

    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    path = mol_list[0].properties.path
    paralel = JobRunner(parallel=True, maxjobs=os.cpu_count())

    # Run MOPAC
    init(path=path, folder='mopac')
    mopac_jobs = list(AMSJob(molecule=mol,
                             settings=mopac_settings(mol.properties.charge),
                             name='MOPAC_'+str(i)) for i, mol in enumerate(mol_list))
    results_dict1 = dict((mopac_job.name, mopac_job.run(paralel)) for mopac_job in mopac_jobs)
    finish()

    # Run COSMO-RS
    init(path=path, folder='crs')
    crs = crs_settings()
    crs_jobs = list(CRSJob(settings=crs, name='COSMO-RS_'+str(i)) for i, mol in enumerate(mol_list))
    for crs_job, results1 in zip(crs_jobs, results_dict1):
        if 'mopac.coskf' in results_dict1[results1].files:
            crs_job.settings.input.Compound._h = '"' + results_dict1[results1]['mopac.coskf'] + '"'
        else:
            crs_job.name = False
    results_dict2 = dict((crs_job.name, crs_job.run(paralel)) for crs_job in crs_jobs if
                         crs_job.name)
    finish()

    # Extract results from the calculations
    Angstrom = Units.convert(1.0, 'Bohr', 'Angstrom')
    for i, mol in enumerate(mol_list):
        prop = mol.properties
        results1 = results_dict1.get('MOPAC_'+str(i))
        results2 = results_dict2.get('COSMO-RS_'+str(i))
        if results1 and 'mopac.coskf' in results1.files:
            prop.surface = KFFile(results1['mopac.coskf']).read('COSMO', 'Area') * Angstrom**2
            prop.volume = KFFile(results1['mopac.coskf']).read('COSMO', 'Volume') * Angstrom**3
            if results2:
                prop.logp = KFFile(results2['$JN.crskf']).read('LOGP', 'logp')[0]

    shutil.rmtree(mopac_jobs[0].path.rsplit('/', 1)[0])
    shutil.rmtree(crs_jobs[0].path.rsplit('/', 1)[0])

    return mol_list


def ams_job_uff_opt(plams_mol, maxiter=2000):
    """
    Runs an AMS UFF constrained geometry optimization.

    plams_mol <plams.Molecule>: The input PLAMS molecule with the 'path', 'name' and
        'mol_indices' properties.
    bonds <list>[<list>[<int>, <float>]]: A nested list of atomic indices and bond orders.
    maxiter <int>: The maximum number of iterations during the geometry optimization.

    return <plams.Molecule>: A PLAMS molecule.
    """
    name = plams_mol.properties.name + '.opt'
    path = plams_mol.properties.path
    mol_indices = plams_mol.properties.indices

    s1 = uff_settings(plams_mol, mol_indices)

    # Run the job (pre-optimization)
    init(path=path, folder=name)
    job = AMSJob(molecule=plams_mol, settings=s1, name=name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    # Delete the PLAMS directory and update MaxIterations
    shutil.rmtree(os.path.join(path, name))
    s1.input.ams.GeometryOptimization.MaxIterations = maxiter

    # Update the atomic coordinates of plams_mol
    for i, atom in enumerate(plams_mol):
        atom.move_to(output_mol[i + 1])

    # Set all H-C=C angles to 120.0 degrees and run the job again
    init(path=path, folder=name)
    job = AMSJob(molecule=fix_h(plams_mol), settings=s1, name=name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], os.path.join(path, name + '.ams.rkf'))
    shutil.copy2(results['uff.rkf'], os.path.join(path, name + '.uff.rkf'))
    shutil.copy2(results[name + '.out'], os.path.join(path, name + '.out'))
    shutil.rmtree(job.path.rsplit('/', 1)[0])

    # Update the atomic coordinates of plams_mol
    for i, atom in enumerate(plams_mol):
        atom.move_to(output_mol[i + 1])

    # Write the reuslts to an .xyz and .pdb file
    plams_mol.properties.name += '.opt'
    export_mol(plams_mol, message='Optimized core + ligands:\t\t')
    plams_mol.properties.name = plams_mol.properties.name.split('.opt')[0]

    return plams_mol
