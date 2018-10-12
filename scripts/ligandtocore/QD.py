import copy
import itertools
import time
import sys

from scm.plams import (Atom, MoleculeError)
from qmflows import molkit

import QD_functions as qd_scripts
import QD_database as qd_database

def prep(arg):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.
    """
    # The start
    time_start = time.time()
    print('\n')

    # Create the result directories
    core_folder, ligand_folder, qd_folder = [qd_scripts.create_dir(name, path=arg['path']) for
                                             name in arg['dir_name_list']]

    # Imports the cores and ligands
    core_dict = qd_scripts.read_mol(arg['input_cores'], core_folder, arg['column'], arg['row'])
    ligand_dict = qd_scripts.read_mol(arg['input_ligands'], ligand_folder, arg['column'],
                                      arg['row'])
   
    # Return the indices of the core dummy atoms
    core_indices = list(prep_core(core_dict[core], core_folder, arg['core_indices'], arg['dummy'], 
                                            arg['core_opt']) for core in core_dict)

    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        database = qd_database.read(ligand_folder, arg['database_name'])
    else:
        database = False

    # Rotate all the ligands and format the resulting list
    ligand_dict = list(prep_ligand(ligand_dict[ligand], ligand, ligand_folder, database, 
                                   arg['ligand_opt'], arg['split']) for ligand in ligand_dict)
    ligand_dict, ligand_indices, database_entries = zip(*ligand_dict)
    ligand_dict = qd_scripts.concatenate_dict(ligand_dict)
    ligand_indices = list(itertools.chain(*ligand_indices))

    # Write new entries to the ligand database
    if arg['use_database']:
        qd_database.write(database_entries, database, ligand_folder)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_dict = list(prep_qd(core_dict[core], core, ligand_dict[ligand], ligand, core_indices[i], 
                           ligand_indices[j], qd_folder) for i, core in enumerate(core_dict) for
                           j, ligand in enumerate(ligand_dict))
    qd_dict, qd_indices = zip(*qd_dict)
    qd_dict = qd_scripts.concatenate_dict(qd_dict)

    # Check if the ADF environment variables are set and optimize the qd with the core frozen
    if arg['qd_opt']:
        sys_var = qd_scripts.check_sys_var()
        if sys_var:
            for i, qd in enumerate(qd_dict):
                qd_scripts.ams_job(qd_dict[qd], qd, qd_folder, qd_indices[i], arg['maxiter'])

    # The End
    time_end = time.time()
    print('\nTotal elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')


def prep_core(core, core_folder, core_indices=[], dummy=0, opt=False):
    """
    Function that handles all core operations.
    """
    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum

    # Optimize the core with RDKit UFF if opt = True. Returns a RDKit molecule
    if opt:
        core = molkit.global_minimum(core)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    if not core_indices:
        core_indices = [(i + 1) for i, atom in enumerate(core.atoms) if atom.atnum == dummy]
        core_indices.reverse()
    else:
        core_indices.sort(reverse=True)
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Set a number of atomic properties
    qd_scripts.set_pdb(core, 'COR', is_core=True)

    # Returns an error if no dummy atoms were found
    if not core_indices:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')
    else:
        return core_indices


def prep_ligand(ligand, ligand_name, ligand_folder, database, opt=True, split=True):
    """
    Function that handles all ligand operations,
    """
    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand, database_entry = qd_scripts.manage_ligand(ligand, ligand_name, ligand_folder, opt, database)

    # Identify functional groups within the ligand and add a dummy atom to the center of mass.
    ligand_dict, ligand_indices = qd_scripts.find_substructure(ligand, ligand_name, split)
    for ligand in ligand_dict:
        ligand_dict[ligand].add_atom(Atom(atnum=0, coords=ligand_dict[ligand].get_center_of_mass()))

    return ligand_dict, ligand_indices, database_entry


def prep_qd(core, core_name, ligand, ligand_name, core_indices, ligand_index, qd_folder):
    """
    Function that handles all quantum dot (qd, i.e. core + all ligands) operations.
    """
    # Rotate and translate all ligands to their position on the core.
    # Returns a list of PLAMS molecules and atomic indices.
    core = copy.deepcopy(core)
    ligand_list = [qd_scripts.rotate_ligand(core, ligand, core_index, ligand_index, i)
                   for i, core_index in enumerate(core_indices)]

    ligand_list, ligand_indices = zip(*ligand_list)
    core.delete_atom(core[-1])

    # Prepare the .pdb filename as a string.
    qd_name = 'Core_' + core_name + '__' + ligand_name

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    qd = qd_scripts.combine_qd(core, ligand_list)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in ligand_indices]
    qd_indices += [i + 1 for i, atom in enumerate(core)]

    qd_scripts.export_mol(qd, qd_folder, qd_name, message='core + ligands:\t\t\t')

    return {qd_name: qd}, qd_indices




argument_dict = {
    # Mandatory arguments: these will have to be manually specified by the user
    'input_cores': {'Cd16Se13.xyz': 'xyz'},
    'input_ligands': {'input_ligands.txt': 'txt'},
    'path': r'D:\QMFlows_DATA',

    # Optional arguments: these can be left to their default values
    'dir_name_list': ['core', 'ligand', 'QD'],
    'smiles_extension': 'txt',
    'column': 0,
    'row': 0,
    'dummy': 'Cl',
    'core_indices': [],
    'ligand_indices': [],
    'database_name': 'ligand_database.xlsx',
    'use_database': True,
    'core_opt': False,
    'ligand_opt': True,
    'qd_opt': False,
    'maxiter': 10000,
    'split': True
}

# For running the script directly from the console
# e.g. python /path_to_QD/QD.py keyword_1:arg_1 keyword_2:arg2 keyword_3:arg3
# SMILES strings should be encapsulated by quotation marks.
# Multiple SMILES string can be passed when seperated by commas, e.g. input_ligands:'OC, OOC, OCCC'
argv = [item for item in sys.argv[1:]]
if argv:
    argv = [item.split(':') for item in argv]
    for item in argv:
        keyword = item[0]
        argument = item[1]
        argument = argument.replace('"', '')
        argument = argument.replace("'", '')
        argument_dict[keyword] = argument.split(',')

# Runs the script: add ligand to core and optimize (UFF) the resulting qd with the core frozen
prep(argument_dict)

"""
input_cores =       The input core(s) as either .xyz, .pdb, .mol, SMILES string, plain text file
                    with SMILES string or a list containing any of the above objects.
input_ligands =     Same as input_cores, except for the ligand(s).
path =              The path where the input and output directories will be saved. Set to
                    os.getcwd() to use the current directory.
dir_name_list =     Names of the to be created directories in path.
smiles_extension =  Extension (without period) of a SMILES string containg plain text file. Relevant
                    if such a file is chosen for input_cores or input_ligands.
column =            The column containing the SMILES strings in the plain text file.
row =               The amount of rows to be ignored in the SMILES string containing column.
                    Should be used when e.g. the first row does not contain a SMILES string
dummy =             The atomic number of atomic symbol of the atoms in the core that should be
                    should be replaced with ligands.
core_indices =      Manually specify the indices of the core dummy atoms instead of utilizing the 
                    'dummy' argument.
ligand_indices =    Manually specifiy the indices of ligand dummy atoms instead of utilizing the
                    find_substructure() function.
database_name =     Name plus extension of the (to be) created ligand database.
use_database =      Export/import results from the (to be) created ligand database. No database will
                    be used and/or maintained when set to False.
core_opt =          Attempt to find the core global minimum using RDKit UFF.
                    WARNING: enabling this will probably ruin the core if care is not taken!
                    Should work fine for organic cores.
ligand_opt =        Attempt to find the ligand global minimum using RDKit UFF.
qd_opt =            Optimize the quantum dot (qd)(i.e core + all ligands) using ADF UFF.
maxiter =           The maximum number of geometry iteration during qd_opt.
split =             Should the ligand be attached to the core in its entirety or should a
                    hydrogen atom/counterion first be removed? Examples are provided below:
                    True:  HO2CR -> -O2CR,  X-.NH4+ -> NH4+  &  Na+.-O2CR -> -O2CR
                    False: HO2CR -> HO2CR,  NH4+ -> NH4+  & -O2CR -> -O2CR
"""
