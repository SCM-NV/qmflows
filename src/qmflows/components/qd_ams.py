__all__ = ['check_sys_var', 'ams_job']

import os
import shutil

from scm.plams import (Settings, AMSJob, init, finish)
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)
# import scm.plams.interfaces.molecule.rdkit as molkit

from .qd_import_export import export_mol
from .qd_functions import adf_connectivity


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
        raise ImportError('No ADF version 2018 detected, aborting ADF job.')


class CRSJob(SCMJob):
    """
    A class for running COSMO-RS jobs.
    """
    _command = 'crs'


def ams_job(plams_mol, maxiter=1000, job='qd_opt'):
    """
    Converts PLAMS connectivity into adf .run script connectivity.

    plams_mol <plams.Molecule>: The input PLAMS molecule.
    maxiter <int>: The maximum number of iterations during the geometry optimization.
    job <str>: The to be run AMS job ('qd_opt', 'qd_sp' or 'ligand_sp')

    return <plams.Molecule>: A PLAMS molecule.
    """
    # Launch an AMS UFF geometry optimization
    if job == 'qd_opt':
        ams_job_uff_opt(plams_mol, maxiter)

    # Launch an MOPAC + COSMO-RS single point
    elif job == 'ligand_sp':
        ams_job_mopac_sp(plams_mol)

    return plams_mol


def ams_job_mopac_sp(plams_mol):
    """
    Runs a MOPAC + COSMO-RS single point.

    plams_mol <plams.Molecule>: The input PLAMS molecule.

    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    source_folder = plams_mol.properties.source_folder

    # MOPAC settings
    mol_name1 = plams_mol.properties.name + '.MOPAC'
    s1 = Settings()
    s1.input.ams.Task = 'SinglePoint'

    # COSMO-RS settings
    mol_name2 = plams_mol.properties.name + '.COSMO-RS'
    s2 = Settings()
    s2.ignore_molecule = True

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

    s2.input.Compound._h = os.path.join(source_folder, mol_name1 + '.t21')

    s2.input.compound._h = os.path.join(os.environ['ADFRESOURCES'], 'ADFCRS/1-Octanol.coskf')
    s2.input.compound.frac1 = 0.725
    s2.input.compound.pvap = 1.01325
    s2.input.compound.tvap = 468.0
    s2.input.compound.meltingpoint = 257.0
    s2.input.compound.flashpoint = 354.0
    s2.input.compound.density = 0.824

    s2.input.compounD._h = os.path.join(os.environ['ADFRESOURCES'], 'ADFCRS/Water.coskf')
    s2.input.compounD.frac1 = 0.275
    s2.input.compounD.frac2 = 1.0
    s2.input.compounD.pvap = 1.01325
    s2.input.compounD.tvap = 373.15
    s2.input.compounD.meltingpoint = 273.15
    s2.input.compounD.hfusion = 1.436
    s2.input.compounD.density = 0.997

    # Run MOPAC
    init(path=source_folder, folder=mol_name1)
    job1 = AMSJob(molecule=plams_mol, settings=s1, name=mol_name1)
    results1 = job1.run()
    finish()

    # Delete the MOPAC directory
    shutil.copy2(results1['ams.rkf'], os.path.join(source_folder, mol_name1 + '.rkf'))
    shutil.rmtree(os.path.join(source_folder, mol_name1))

    # Run COSMO-RS
    init(path=source_folder, folder=mol_name2)
    job2 = AMSJob(settings=s2, name=mol_name2)
    results2 = job2.run()
    finish()

    # Read logp
    with open(results2['AcOH.out']) as file:
        file = file.read().splitlines()
    logp_index = file.index(' Infinite dilute Partition coefficient') + 2
    logp = file[logp_index].split()[-1]

    # Delete the PLAMS directory
    shutil.copy2(results2[mol_name2 + '.out'], os.path.join(source_folder, mol_name2 + '.out'))
    shutil.rmtree(os.path.join(source_folder, mol_name2))

    # Create three new properties for plams_mol
    plams_mol.properties.surface = 0.0
    plams_mol.properties.volume = 0.0
    plams_mol.properties.logp = logp

    return plams_mol


def ams_job_uff_opt(plams_mol, maxiter=2000):
    """
    Runs an AMS UFF constrained geometry optimization.

    plams_mol <plams.Molecule>: The input PLAMS molecule.
    bonds <list>[<list>[<int>, <float>]]: A nested list of atomic indices and bond orders.
    maxiter <int>: The maximum number of iterations during the geometry optimization.

    return <plams.Molecule>: A PLAMS molecule.
    """
    mol_name = plams_mol.properties.name + '.opt'
    source_folder = plams_mol.properties.source_folder
    mol_indices = plams_mol.properties.qd_indices

    # AMS settings (UFF constrained geometry optimization)
    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Constraints.Atom = mol_indices
    s.input.ams.System.BondOrders._1 = adf_connectivity(plams_mol)
    s.input.ams.GeometryOptimization.MaxIterations = maxiter
    s.input.uff.Library = 'UFF'

    # Run the job
    init(path=source_folder, folder=mol_name)
    job = AMSJob(molecule=plams_mol, settings=s, name=mol_name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    """
    # Run the job
    init(path=source_folder, folder=mol_name)
    job = AMSJob(molecule=output_mol, settings=s, name=mol_name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()
    """
    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], os.path.join(source_folder, mol_name + '.ams.rkf'))
    shutil.copy2(results['uff.rkf'], os.path.join(source_folder, mol_name + '.uff.rkf'))
    shutil.copy2(results[mol_name + '.out'], os.path.join(source_folder, mol_name + '.out'))
    #shutil.rmtree(os.path.join(source_folder, mol_name))

    # Update the atomic coordinates of plams_mol
    for i, atom in enumerate(plams_mol):
        atom.move_to(output_mol[i + 1])

    # Write the reuslts to an .xyz and .pdb file
    plams_mol.properties.name += '.opt'
    export_mol(plams_mol, message='Optimized core + ligands:\t\t')

    return output_mol
