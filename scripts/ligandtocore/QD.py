import copy
import itertools
import time
import sys

from scm.plams import (Atom, MoleculeError, Settings)
from qmflows import molkit

import QD_functions as QD_scripts
import QD_database
import QD_import_export as QD_inout


def prep(input_ligands, input_cores, path, arg):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.
    """
    # The start
    time_start = time.time()
    print('\n')

    # Create the result directories (if they do not exist), ligands and cores
    folder_list = [QD_inout.create_dir(name, path) for name in arg['dir_name_list']]
    ligand_list = QD_inout.read_mol(input_ligands, folder_list[1], arg['column'], arg['row'])
    core_list = QD_inout.read_mol(input_cores, folder_list[0], arg['column'], arg['row'],
                                  is_core=True)

    # Adds the indices of the core dummy atoms to core.properties.core
    for core in core_list:
        prep_core(core, arg['core_indices'], arg['dummy'], arg['core_opt'])

    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        database = QD_database.read(folder_list[1], arg['database_name'])
    else:
        database = False

    # Optimize all ligands and find their functional groups
    ligand_list = list(prep_ligand(ligand, database, arg['ligand_opt'], arg['split']) for
                       ligand in ligand_list)
    ligand_list = list(itertools.chain(*ligand_list))

    # Write new entries to the ligand database
    if arg['use_database']:
        QD_database.write(ligand_list, database)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_list = list(prep_qd(core, ligand, folder_list[2]) for core in core_list for
                   ligand in ligand_list)

    # Check if the ADF environment variables are set and optimize the qd with the core frozen
    if arg['qd_opt']:
        if QD_scripts.check_sys_var():
            for qd in qd_list:
                QD_scripts.ams_job(qd, arg['maxiter'])

    # The End
    time_end = time.time()
    print('\nTotal elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')

    return qd_list


def prep_core(core, core_indices, dummy=0, opt=False):
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
    core_indices.sort(reverse=True)
    core.properties.core_indices = core_indices
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Returns an error if no dummy atoms were found
    if not core_indices:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')


def prep_ligand(ligand, database, opt=True, split=True):
    """
    Function that handles all ligand operations,
    """
    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand = QD_scripts.optimize_ligand(ligand, opt, database)

    # Identify functional groups within the ligand and add a dummy atom to the center of mass.
    ligand_list = QD_scripts.find_substructure(ligand, split)

    return ligand_list


def prep_qd(core, ligand, qd_folder):
    """
    Function that handles all quantum dot (qd, i.e. core + all ligands) operations.
    """
    # Rotate and translate all ligands to their position on the core.
    # Returns a list of PLAMS molecules and atomic indices.
    core = copy.deepcopy(core)
    ligand_list = [QD_scripts.rotate_ligand(core, ligand, i, index) for i, index in
                   enumerate(core.properties.core_indices)]
    ligand_list, ligand_indices = zip(*ligand_list)
    core.delete_atom(core[-1])

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    qd = QD_scripts.combine_qd(core, ligand_list)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in ligand_indices]
    qd_indices += [i + 1 for i, atom in enumerate(core)]

    qd.properties = Settings()
    qd.properties.qd_indices = qd_indices
    qd.properties.name = core.properties.name + '__' + ligand.properties.name
    qd.properties.source_folder = qd_folder
    QD_inout.export_mol(qd, message='core + ligands:\t\t\t')

    return qd


# Mandatory arguments: these will have to be manually specified by the user
# Key: filename
# Argument: string containing the filetype (i.e. 'xyz', 'pdb', 'mol', 'smiles', 'folder', 'txt')
# Or argument: list containing [0] the filetype (see above) and [1] if bonds should be guessed used (Bool)
# By default guess_bonds() is only enabled for .xyz files

input_cores = {
        'Cd16Se13.xyz': ['xyz', False]
        }

input_ligands = {
        'input_ligands.txt': 'txt'
        }

path = r'D:\QMFlows_DATA'

# Optional arguments: these can be left to their default values
argument_dict = {
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
qd_list = prep(input_ligands, input_cores, path, argument_dict)

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
