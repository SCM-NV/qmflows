__all__ = ['check_sys_var', 'ams_job_mopac_sp', 'qd_opt']

import math
import os
import shutil

import numpy as np

from scm.plams import (Settings, AMSJob, init, finish, Units, add_to_class, AMSResults)
from scm.plams.tools.kftools import KFFile
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)
from scm.plams.core.jobrunner import JobRunner
import scm.plams.interfaces.molecule.rdkit as molkit

from .qd_import_export import export_mol
from .qd_functions import (adf_connectivity, fix_h, fix_carboxyl, get_time)


def check_sys_var():
    """
    Check if all ADF environment variables are set and if the 2018 version of ADF is installed.
    """
    sys_var = ['ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            print(get_time() +
                  'WARNING: The environment variable ' + sys_var[i] + ' has not been set')
    if False in sys_var_exists:
        raise EnvironmentError(get_time() + 'One or more ADF environment variables have '
                               'not been set, aborting ADF job.')
    if '2018' not in os.environ['ADFHOME']:
        error = get_time() + 'No ADF version 2018 detected in ' + os.environ['ADFHOME']
        error += ', aborting ADF job.'
        raise ImportError(error)


def qd_opt(mol, database, arg):
    """
    Check if the to be optimized quantom dot has previously been optimized.
    Pull if the structure from the database if it has, otherwise perform a geometry optimization.
    mol <plams.Molecule> The input quantom dot with the 'name' property.
    database <pd.DataFrame>: A database of previous calculations.
    return <plams.Molecule>: An optimized quantom dot.
    """
    name = mol.properties.name.rsplit('.', 1)[0]
    if database is None:
        mol = ams_job_uff_opt(mol, arg['maxiter'])
        mol.properties.entry = False
    elif database.empty or name not in list(database['Quantum_dot_name']):
        mol = ams_job_uff_opt(mol, arg['maxiter'])
        mol.properties.entry = True
    else:
        index = list(database['Quantum_dot_name']).index(name)
        try:
            mol_new = molkit.readpdb(database['Quantum_dot_opt_pdb'][index])
            mol_new.properties = mol.properties
            mol = mol_new
        except FileNotFoundError:
            mol = ams_job_uff_opt(mol, arg['maxiter'])
            mol.properties.entry = True

    return mol


class CRSResults(SCMResults):
    """
    A class for accessing results of COSMO-RS jobs.
    """
    _kfext = '.crskf'
    _rename_map = {'CRSKF': '$JN.crskf'}

    def get_solvation_energy(self, kf):
        """
        Returns the solvation energy from a Activity Coefficients calculation.
        """
        return KFFile(self[kf]).read('ACTIVITYCOEF', 'deltag')[0]


class CRSJob(SCMJob):
    """
    A class for running COSMO-RS jobs.
    """
    _command = 'crs'
    _result_type = CRSResults


@add_to_class(AMSResults)
def get_entropy(self, freqs, T=298.15):
    """
    Calculate the translational, vibrational and rotational entropy.
    All units and constants are in SI units.
    self <plams.AMSResults>: An AMSResults object of a vibrational analysis
    freqs <np.ndarray>: A numpy array containg vibrational frequencies in atomic units
    T <float>: The temperature in Kelvin
    Return <np.ndarray>: A numpy array containing the translational, rotational and
        vibrational contributions to the entropy
    """
    # Define constants
    kT = 1.380648 * 10**-23 * T  # Boltzmann constant * temperature
    h = 6.6260701 * 10**-34  # Planck constant
    hv_kT = (h * freqs) / kT  # (Planck constant * frequencies) / (Boltzmann * temperature)
    R = 8.31445  # Gas constant
    V_Na = ((R * T) / 10**5) / Units.constants['NA']  # Volume(1 mol ideal gas) / Avogadro's number
    pi = math.pi

    # Extract atomic masses and coordinates
    mol = self.get_main_molecule()
    m = np.array([at.mass for at in mol]) * 1.6605390 * 10**-27
    x, y, z = mol.to_array().T * 10**-10

    # Calculate the rotational partition function
    inertia = np.array([sum(m*(y**2 + z**2)), -sum(m*x*y), -sum(m*x*z),
                        -sum(m*x*y), sum(m*(x**2 + z**2)), -sum(m*y*z),
                        -sum(m*x*z), -sum(m*y*z), sum(m*(x**2 + y**2))]).reshape(3, 3)
    inertia = np.product(np.linalg.eig(inertia)[0])
    q_rot = pi**0.5 * ((8 * pi**2 * kT) / h**2)**1.5 * inertia**0.5

    # Calculate the translational, rotational and vibrational entropy (divided by R)
    S_trans = 1.5 + np.log(V_Na * ((2 * pi * sum(m) * kT) / h**2)**1.5)
    S_rot = 1.5 + np.log(q_rot)
    S_vib = sum(hv_kT / np.expm1(hv_kT) - np.log(1 - np.exp(-hv_kT)))

    return R * np.array([S_trans, S_rot, S_vib])


@add_to_class(AMSResults)
def get_thermo(self, kf, T=298.15, export=['E', 'G'], unit='kcal/mol'):
    """
    Extract and return Gibbs free energies, entropies and/or enthalpies from an AMS KF file.
    All vibrational frequencies smaller than 100 cm**-1 are set to 100 cm**-1.
    self <plams.AMSResults>: An AMSResults object resulting from a vibrational analysis
    kf <str>: The name+extension of the AMS KF file containing energies and frequencies
    T <float>: The temperature in Kelvin
    export <list>[<str>]: An iterable containing strings of the to be exported energies:
        'E': Electronic energy
        'U': Interal energy (E + U_nuc)
        'H': Enthalpy (U + pV)
        'S': Entropy
        'G': Gibbs free energy (H - T*S)
    unit <str>: The unit of the to be returned energies.
    Return <float> or <dict>[<float>]: An energy or dictionary of energies
    """
    # Get frequencies; set all frequencies smaller than 100 cm**-1 to 100 cm**-1
    freqs = np.array(KFFile(self[kf]).read('Vibrations', 'Frequencies[cm-1]'))
    freqs[freqs < 100] = 100
    freqs *= 100 * Units.constants['c']

    # hv_kT = (Planck constant * frequencies) / (Boltzmann constant * temperature)
    hv_kT = (6.6260701 * 10**-34 * freqs) / (1.380648 * 10**-23 * T)
    RT = 8.31445 * T  # Gas constant * temperature

    # Extract and/or calculate the various energies
    E = KFFile(self[kf]).read('AMSResults', 'Energy') * Units.conversion_ratio('Hartree', 'kj/mol')
    E *= 1000
    U = E + RT * (3.0 + sum(0.5 * hv_kT + hv_kT / np.expm1(hv_kT)))
    H = U + RT
    S = sum(self.get_entropy(freqs, T=T))
    G = H - T * S

    ret = {'E': E, 'U': U, 'H': H, 'S': S, 'G': G}

    if len(export) == 1:
        return Units.convert(ret[export[0]], 'kj/mol', unit) / 1000
    return {i: Units.convert(ret[i], 'kj/mol', unit) / 1000 for i in ret if i in export}


def uff_settings(plams_mol, mol_indices, maxiter=2000):
    """
    UFF settings for a constrained geometry optimization.
    """
    s1 = Settings()
    s1.pickle = False

    s1.input.ams.Task = 'GeometryOptimization'
    s1.input.ams.Constraints.Atom = mol_indices
    s1.input.ams.System.BondOrders._1 = adf_connectivity(plams_mol)
    s1.input.ams.GeometryOptimization.MaxIterations = int(maxiter / 2)
    s1.input.uff.Library = 'UFF'

    return s1


def crs_settings(solute, solvent):
    """
    COSMO-RS settings for Activity Coefficient calculations using 'MOPAC PM6' parameters.
    Yields the solvation energy, among other things.
    """
    s2 = Settings()
    s2.ignore_molecule = True
    s2.pickle = False

    s2.input.Property._h = 'activitycoef'

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

    s2.input.Compound._h = '"' + solute + '"'

    path = os.path.join(os.path.dirname(__file__), 'coskf')
    s2.input.compound._h = '"' + os.path.join(path, solvent + '.coskf') + '"'
    s2.input.compound.frac1 = 1.0

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


def ams_job_mopac_sp(mol):
    """
    Runs a MOPAC + COSMO-RS single point.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    path = mol.properties.path
    angstrom = Units.convert(1.0, 'Bohr', 'Angstrom')
    solvents = ('Acetone', 'Acetonitrile', 'DMF', 'DMSO', 'Ethanol',
                'EtOAc', 'Hexane', 'Toluene', 'Water')

    # Run MOPAC
    init(path=path, folder='ligand')
    mopac_job = AMSJob(molecule=mol,
                       settings=mopac_settings(mol.properties.charge),
                       name='MOPAC')
    mopac_results = mopac_job.run()

    if 'mopac.coskf' in mopac_results.files:
        # Extract properties from mopac.coskf
        coskf = mopac_results['mopac.coskf']
        mol.properties.surface = KFFile(coskf).read('COSMO', 'Area') * angstrom**2
        mol.properties.volume = KFFile(coskf).read('COSMO', 'Volume') * angstrom**3

        # Run COSMO-RS in parallel
        parallel = JobRunner(parallel=True)
        crs_jobs = [CRSJob(settings=crs_settings(coskf, solv), name='COSMO-RS_'+solv) for
                    solv in solvents]
        crs_result = [job.run(jobrunner=parallel) for job in crs_jobs]

        # Extract properties from COSMO-RS_solv.crskf
        solvation_energy = Settings()
        for results, solv in zip(crs_result, solvents):
            if 'COSMO-RS_' + solv + '.crskf' in results.files:
                solvation_energy[solv] = results.get_solvation_energy('$JN.crskf')
            else:
                solvation_energy[solv] = None
        mol.properties.solvation = solvation_energy
    finish()

    shutil.rmtree(mopac_job.path.rsplit('/', 1)[0])

    return mol


def ams_job_uff_opt(mol, maxiter=2000):
    """
    Runs an AMS UFF constrained geometry optimization.
    mol <plams.Molecule>: The input PLAMS molecule with the 'path', 'name' and
        'mol_indices' properties.
    bonds <list>[<list>[<int>, <float>]]: A nested list of atomic indices and bond orders.
    maxiter <int>: The maximum number of iterations during the geometry optimization.
    return <plams.Molecule>: A PLAMS molecule.
    """
    name = mol.properties.name + '.opt'
    path = mol.properties.path
    mol_indices = mol.properties.indices

    # Run the job (pre-optimization)
    s1 = uff_settings(mol, mol_indices, maxiter=maxiter)
    init(path=path, folder='Quantum_dot')
    job = AMSJob(molecule=mol, settings=s1, name='UFF_part1')
    results = job.run()
    mol.from_iterable(results.get_main_molecule())

    # Fix all O=C-O and H-C=C angles and continue the job
    s1.input.ams.Properties.NormalModes = 'Yes'
    job = AMSJob(molecule=fix_carboxyl(fix_h(mol)), settings=s1, name='UFF_part2')
    results = job.run()
    mol.from_iterable(results.get_main_molecule())
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], os.path.join(path, name + '.ams.rkf'))
    shutil.copy2(results['uff.rkf'], os.path.join(path, name + '.uff.rkf'))
    shutil.copy2(results['$JN.out'], os.path.join(path, name + '.out'))
    shutil.rmtree(job.path.rsplit('/', 1)[0])

    # Write the reuslts to an .xyz and .pdb file
    mol.properties.name += '.opt'
    export_mol(mol, message='Optimized core + ligands:\t\t')
    mol.properties.name = mol.properties.name.split('.opt')[0]

    return mol
