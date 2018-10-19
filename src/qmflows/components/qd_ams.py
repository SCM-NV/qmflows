__all__ = ['check_sys_var', 'ams_job']

import os
import shutil

from scm.plams import (Settings, AMSJob, init, finish)
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit.Chem import Bond

from .qd_import_export import export_mol


def check_sys_var():
    """
    Check if all ADF environment variables are set.
    """
    sys_var = ['ADFBIN', 'ADFHOME', 'ADFRESOURCES', 'SCMLICENSE']
    sys_var_exists = [item in os.environ for item in sys_var]
    for i, item in enumerate(sys_var_exists):
        if not item:
            print('WARNING: The environment variable ' + sys_var[i] + ' has not been set')
    if False in sys_var_exists:
        raise EnvironmentError('One or more ADF environment variables have not been set, aborting '
                               'ADF job.')
    version = float(os.environ['ADFHOME'].rsplit('/').rsplit('\\').split('ADF')[-1])
    if version < 2018:
        raise ImportError('ADF version', version, 'detected.'
                          'ADF2018 or later is required, aborting ADF job.')


def ams_job(plams_mol, maxiter=1000, job='qd_opt'):
    """
    Converts PLAMS connectivity into adf .run script connectivity.
    """
    if job == 'qd_opt' or job == 'qd_int':
        # Create a list of aromatic bond indices
        rdmol = molkit.to_rdmol(plams_mol)
        aromatic = [Bond.GetIsAromatic(bond) for bond in rdmol.GetBonds()]
        aromatic = [i for i, item in enumerate(aromatic) if item]

        # Create a connectivity list; aromatic bonds get a bond order of 1.5
        at1 = [plams_mol.atoms.index(bond.atom1) + 1 for bond in plams_mol.bonds]
        at2 = [plams_mol.atoms.index(bond.atom2) + 1 for bond in plams_mol.bonds]
        bonds = [bond.order for bond in plams_mol.bonds]
        for i, bond in enumerate(plams_mol.bonds):
            if i in aromatic:
                bonds[i] = 1.5
        bonds = [str(at1[i]) + ' ' + str(at2[i]) + ' ' + str(bond) for i, bond in enumerate(bonds)]

        # Launch an AMS UFF constrained geometry optimization
        if job == 'qd_opt':
            plams_mol = ams_job_uff_opt(plams_mol, bonds, maxiter)

        # Launch an AMS UFF single point
        if job == 'qd_sp':
            plams_mol = ams_job_uff_sp(plams_mol, bonds)

    # Launch an MOPAC + COSMO-RS single point
    if job == 'ligand_sp':
        plams_mol = ams_job_mopac_sp(plams_mol, maxiter)

    return plams_mol


def ams_job_mopac_sp(plams_mol, maxiter):
    """
    Runs an MOPAC + COSMO-RS single point.
    """
    source_folder = plams_mol.properties.source_folder

    # MOPAC settings
    mol_name1 = plams_mol.properties.name + '.MOPAC'
    s1 = Settings()
    s1.input.ams.Task = 'SinglePoint'
    s1.input.ams.GeometryOptimization.MaxIterations = maxiter

    # COSMO-RS settings
    mol_name2 = plams_mol.properties.name + '.COSMO-RS'
    s2 = Settings()
    'To be added'

    # Run MOPAC
    init(path=source_folder, folder=mol_name1)
    job1 = AMSJob(molecule=plams_mol, settings=s1, name=mol_name1)
    results1 = job1.run()
    finish()

    # Run COSMO-RS
    init(path=source_folder, folder=mol_name2)
    job2 = AMSJob(molecule=plams_mol, settings=s2, name=mol_name2)
    results2 = job2.run()
    finish()

    # Delete the PLAMS directories
    shutil.rmtree(os.path.join(source_folder, mol_name1))
    shutil.rmtree(os.path.join(source_folder, mol_name2))

    return plams_mol


def ams_job_uff_opt(plams_mol, bonds, maxiter):
    """
    Runs an AMS UFF constrained geometry optimization.
    """
    mol_name = plams_mol.properties.name + '.opt'
    source_folder = plams_mol.properties.source_folder
    mol_indices = plams_mol.properties.qd_indices

    # General AMS settings
    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Constraints.Atom = mol_indices
    s.input.ams.System.BondOrders._1 = bonds
    s.input.ams.GeometryOptimization.MaxIterations = maxiter
    s.input.uff.Library = 'UFF'

    # Run the job
    init(path=source_folder, folder=mol_name)
    job = AMSJob(molecule=plams_mol, settings=s, name=mol_name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], os.path.join(source_folder, mol_name + '.rkf'))
    shutil.copy2(results[mol_name + '.out'], os.path.join(source_folder, mol_name + '.out'))
    shutil.rmtree(os.path.join(source_folder, mol_name))

    # Update the atomic coordinates of plams_mol
    for i, atom in enumerate(plams_mol):
        atom.move_to(output_mol[i + 1])

    # Write the reuslts to an .xyz and .pdb file
    plams_mol.properties.name += '.opt'
    export_mol(plams_mol, message='Optimized core + ligands:\t\t')

    return output_mol


def ams_job_uff_sp(plams_mol, bonds, maxiter):
    """
    Runs an AMS UFF single point.
    """
    mol_name = plams_mol.properties.name
    source_folder = plams_mol.properties.source_folder

    # General AMS settings
    s = Settings()
    s.input.ams.Task = 'SinglePoint'
    s.input.ams.System.BondOrders._1 = bonds
    s.input.uff.Library = 'UFF'

    # Run the job
    init(path=source_folder, folder=mol_name)
    job = AMSJob(molecule=plams_mol, settings=s, name=mol_name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    # Delete the PLAMS directory
    shutil.rmtree(os.path.join(source_folder, mol_name))

    return output_mol

