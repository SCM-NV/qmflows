__all__ = ['check_sys_var', 'ams_job_mopac_crs', 'ams_job_mopac_opt', 'ams_job_mopac_sp',
           'ams_job_uff_opt']

import math
import shutil
import os
from os.path import (dirname, join)

import numpy as np

from scm.plams.core.functions import (init, finish, add_to_class, config)
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.settings import Settings
from scm.plams.tools.kftools import KFFile
from scm.plams.tools.units import Units
from scm.plams.interfaces.adfsuite.ams import (AMSJob, AMSResults)
from scm.plams.interfaces.adfsuite.scmjob import (SCMJob, SCMResults)

from ..templates.templates import get_template
from .qd_import_export import export_mol
from .qd_functions import (adf_connectivity, fix_h, fix_carboxyl, get_time, from_plams_mol)


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
    x, y, z = mol.as_array().T * 10**-10

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


def ams_job_mopac_sp(mol):
    """
    Runs a gas-phase MOPAC single point.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    # Run MOPAC
    config.default_jobmanager.settings.hashing = None
    s = get_template('qd.json')['MOPAC single point']
    job = AMSJob(molecule=mol, settings=s, name='MOPAC')
    results = job.run()
    results.wait()
    mol.properties.energy.E = KFFile(results['mopac.rkf']).read('AMSResults', 'Energy')

    return mol


def ams_job_mopac_opt(mol):
    """
    Runs a gas-phase MOPAC geometry optimization.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    # Run MOPAC
    s = get_template('qd.json')['MOPAC geometry optimization']
    job = AMSJob(molecule=mol, settings=s, name='MOPAC')
    results = job.run()
    results.wait()
    mol.properties.energy.E = KFFile(results['mopac.rkf']).read('AMSResults', 'Energy')

    return mol


def ams_job_mopac_crs(mol):
    """
    Runs a MOPAC + COSMO-RS single point.
    mol <plams.Molecule>: A PLAMS molecule with the 'path' and 'charge' properties.
    return <plams.Molecule>: a PLAMS molecule with the surface, volume and logp properties.
    """
    path = mol.properties.path
    angstrom = Units.conversion_ratio('Bohr', 'Angstrom')
    solvents = ('Acetone', 'Acetonitrile', 'DMF', 'DMSO', 'Ethanol',
                'EtOAc', 'Hexane', 'Toluene', 'Water')
    solv_path = join(join(dirname(dirname(__file__)), 'data'), 'coskf')

    # Run MOPAC
    init(path=path, folder='ligand')
    s1 = get_template('qd.json')['COSMO-MOPAC single point']
    s1.input.ams.System.Charge = mol.properties.charge
    mopac_job = AMSJob(molecule=mol, settings=s1, name='MOPAC')
    mopac_results = mopac_job.run()

    if 'mopac.coskf' in mopac_results.files:
        # Extract properties from mopac.coskf
        coskf = mopac_results['mopac.coskf']
        mol.properties.surface = KFFile(coskf).read('COSMO', 'Area') * angstrom**2
        mol.properties.volume = KFFile(coskf).read('COSMO', 'Volume') * angstrom**3

        # Prepare COSMO-RS a list of settings; one for each solvent
        parallel = JobRunner(parallel=True)
        s2 = get_template('qd.json')['COSMO-RS activity coefficient']
        s2.input = get_template('crs.json')['MOPAC PM6']
        s2.input.Compound._h = '"' + coskf + '"'
        s2_list = []
        for solv in solvents:
            s2_tmp = s2.copy()
            s2_tmp.input.compound._h = '"' + join(solv_path, solv + '.coskf') + '"'
            s2_list.append(s2_tmp)

        # Run COSMO-RS in parallel
        crs_jobs = [CRSJob(settings=s2, name='COSMO-RS_'+solv) for s2 in s2_list]
        crs_result = [job.run(jobrunner=parallel) for job in crs_jobs]

        # Extract properties from COSMO-RS_solv.crskf
        solvation_energy = Settings()
        for results, solv in zip(crs_result, solvents):
            results.wait()
            if 'COSMO-RS_' + solv + '.crskf' in results.files:
                solvation_energy[solv] = results.get_solvation_energy('$JN.crskf')
            else:
                solvation_energy[solv] = None
        mol.properties.energy = solvation_energy
    finish()

    if len(mol.properties.solvation) == len(solvents):
        shutil.rmtree(mopac_job.path.rsplit('/', 1)[0])

    return mol


def ams_job_uff_opt(mol, maxiter=1000, get_freq=False, fix_angle=True):
    """
    Runs an AMS UFF constrained geometry optimization.
    mol <plams.Molecule>: The input PLAMS molecule with the 'path', 'name' and
        'mol_indices' properties.
    bonds <list>[<list>[<int>, <float>]]: A nested list of atomic indices and bond orders.
    maxiter <int>: The maximum number of iterations during the geometry optimization.
    get_freq <bool>: Perform a frequency analyses after the geometry optimization.
    return <plams.Molecule>: A PLAMS molecule.
    """
    name = mol.properties.name + '.opt'
    path = mol.properties.path
    constrains = mol.properties.indices

    # Prepare the job settings
    init(path=path, folder='Quantum_dot')
    s = get_template('qd.json')['UFF constrained optimization']
    s.input.ams.Constraints.Atom = constrains
    s.input.ams.System.BondOrders._1 = adf_connectivity(mol)
    s.input.ams.GeometryOptimization.MaxIterations = int(maxiter / 2)
    if get_freq:
        s.input.ams.Properties.NormalModes = 'Yes'

    # Run the job (pre-optimization)
    job = AMSJob(molecule=mol, settings=s, name='UFF_part1')
    results = job.run()
    output_mol = results.get_main_molecule()
    mol.from_plams_mol(output_mol)
    if get_freq:
        results.wait()
        mol.properties.energy = Settings(results.get_thermo('uff.rkf'))

    # Fix all O=C-O and H-C=C angles and continue the job
    if fix_angle:
        job = AMSJob(molecule=fix_carboxyl(fix_h(mol)), settings=s, name='UFF_part2')
        results = job.run()
        output_mol = results.get_main_molecule()
        mol.from_plams_mol(output_mol)
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], join(path, name + '.ams.rkf'))
    shutil.copy2(results['uff.rkf'], join(path, name + '.uff.rkf'))
    shutil.copy2(results['$JN.out'], join(path, name + '.out'))
    shutil.rmtree(job.path.rsplit('/', 1)[0])

    # Write the reuslts to an .xyz and .pdb file
    mol.properties.name += '.opt'
    export_mol(mol, message='Optimized core + ligands:\t\t')
    mol.properties.name = mol.properties.name.split('.opt')[0]

    return mol
