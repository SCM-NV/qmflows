from scm.plams import *
from qmflows import molkit
import itertools
import time
import copy
import numpy as np
import os
import QD_functions as QD


def prep_core(core, core_folder, dummy=0, opt=True):
    """
    function that handles are core operations
    Identify the ligand anchoring sites on the core and find the center of mass
    """    
    # Checks the if the dummy atom ligand placeholder is provided by its atomic number (int) or atomic symbol (string)
    # Returns an error if neither an integer nor string is provided
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum
    
    if opt:
        core = QD.global_minimum(core, core_folder)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core 
    # An additional dummy atom is added at the core center of mass for the purpose of orietating the ligands in QD.find_substructure(...)
    core_indices = [(i + 1) for i, atom in enumerate(core.atoms) if atom.atnum == dummy]
    core_indices.reverse()
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    if len(core_indices) == 0:
        raise MoleculeError(Atom(atnum=dummy).symbol + ' was specified as dummy atom, yet no dummy atoms were found on the core')
    else:
        return core_indices


def prep_ligand(ligand, ligand_folder, database, opt=True):
    """
    function that handles are ligand operations
    """
    # checks if the database exists and if the ligand is already present
    # returns a list of booleans
    if database:
        matches = [molkit.to_rdmol(ligand).HasSubstructMatch(mol) for mol in database[4]]
    else:
        matches = [False]

    # if the ligand is already present in the database: append the structure, else optimize the ligand and create a new entry for the database
    if any(matches):
        index = matches.index(True)
        ligand = molkit.readpdb(os.path.join(ligand_folder, str(database[3][index])))
        database_entry = False
    else:
        if opt:
            ligand = QD.global_minimum(ligand, ligand_folder)
        database_entry = QD.create_entry(ligand, opt)
        
    # find functional groups
    ligand.add_atom(Atom(atnum=0, coords=(ligand.get_center_of_mass())))
    ligand_list = QD.find_substructure(ligand)
    ligand_list = list(ligand_list)
    ligand_list.append(database_entry)
    
    return ligand_list


def prep_core_ligand(core, ligand, core_indices, ligand_index, core_ligand_folder):
    """
    function that handles all core+ligand operations
    add all ligands to the core
    """
    # Rotate and translate all ligands to their position on the core
    # Returns a list with sublist [0] containing the rotated ligands (PLAMS Molecules) and [1] the heteroatoms (PLAMS Atoms) of the rotated ligands to be attached to the core
    # All core dummy atoms are deleted
    core = copy.deepcopy(core)
    ligand_list = [QD.rotate_ligand(core, ligand, core_index, ligand_index) for core_index in core_indices]
    ligand_list = np.array(ligand_list).T.tolist()

    core.delete_atom(core[-1])

    # Prepare the .pdb filename (string)
    core_name = core.get_formula()
    ligand_name = ligand_list[0][0].get_formula() + '_@_' + ligand[ligand_index + 1].symbol + str(ligand_index + 1)
    pdb_name = str('core_' + core_name + '__and__ligand_' + ligand_name)

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule)
    core_ligand = QD.combine_core_ligand(core, ligand_list[0])
    
    # Create list with sublists containg [0] ligand heteroatom indices (integer), [1] core selenium indices (integer) and [2] core cadmium indices (integer)
    atnum_list = []
    for atom in core:
        if atom.atnum not in atnum_list:
            atnum_list.append(atom.atnum)

    core_ligand_indices = [[core_ligand.atoms.index(atom) + 1 for atom in ligand_list[1]]]
    core_ligand_indices += [[i + 1 for i,atom in enumerate(core) if atom.atnum == atnum] for atnum in atnum_list]
    core_ligand_indices = [item for sublist in core_ligand_indices for item in sublist]
    
    # Assign residue names and formal atomic charges
    #core_ligand = molkit.to_rdmol(core_ligand)
    #core_ligand = QD.prepare_pdb(core_ligand, core, ligand, core_ligand_indices)
    
    # Optimize the combined core and ligands, writing the results to a .pdb file
    molkit.writepdb(core_ligand, os.path.join(core_ligand_folder, pdb_name + '.pdb'))
    print('core + ligands:\t\t\t' + pdb_name + '.pdb')
    #core_ligand = molkit.from_rdmol(core_ligand)
    #molkit.writepdb(core_ligand, os.path.join(core_ligand_folder, pdb_name + '.opt.pdb'))
    #print('\nOptimized core + ligands:\t' + pdb_name + '.opt.pdb')

    return core_ligand, pdb_name, core_ligand_indices





# The start
time_start = time.time()
print('\n')

# managing the result directories
dir_name_list = ['core', 'ligand', 'core_ligand']
dir_path_list = [QD.create_dir(name, path=r'/Users/basvanbeek/Documents/CdSe/Week_5') for name in dir_name_list]
core_folder = dir_path_list[0]
ligand_folder = dir_path_list[1]
core_ligand_folder = dir_path_list[2]

# Accepted inputs: .xyz/.pdb file, SMILES string, plain text file with SMILES strings or a list of aforementioned objects
input_cores = ['Cd176Se147_Cl58.xyz']
input_ligands = ['OC(CCC)C(CCCCCCC)CCCCCCC']

# Imports the cores and ligands
core_list = QD.read_mol(dir_path_list[0], input_cores)
ligand_list = QD.read_mol(dir_path_list[1], input_ligands, smiles_column=0, smiles_extension='.txt')

# Return the indices of the core dummy atoms
core_indices = [prep_core(core, core_folder, dummy='Cl', opt=False) for core in core_list]

# Open the ligand database and check if the specified ligand(s) is already present
database = QD.read_database(ligand_folder)
ligand_list = [prep_ligand(ligand, ligand_folder, database, opt=True) for ligand in ligand_list]

# formating of ligand_list
database_entries = [item[2] for item in ligand_list]
ligand_indices = [item[1] for item in ligand_list]
ligand_indices = list(itertools.chain(*ligand_indices))
ligand_list = [item[0] for item in ligand_list]
ligand_list = list(itertools.chain(*ligand_list))

# writing new entries to the ligand database
QD.write_database(database_entries, ligand_folder, database)

# combine the core with the ligands, yielding core_ligand
core_ligand_list = [prep_core_ligand(core, ligand, core_indices[i], ligand_indices[j], core_ligand_folder) for i,core in enumerate(core_list) for j,ligand in enumerate(ligand_list)]

# formating of core_ligand_list
core_ligand_indices = [item[2] for item in core_ligand_list]
pdb_name_list = [item[1] for item in core_ligand_list]
core_ligand_list = [item[0] for item in core_ligand_list]

# optimize core_ligand with the core frozen
[QD.run_ams_job(core_ligand, pdb_name_list[i], core_ligand_folder, core_ligand_indices[i]) for i,core_ligand in enumerate(core_ligand_list)]

# The End
time_end = time.time()
print("\nTotal elapsed time:\t\t" + "%.4f" % (time_end - time_start) + ' sec')
