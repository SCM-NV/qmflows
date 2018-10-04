from scm.plams import *
from qmflows import molkit
import itertools
import time
import copy
import numpy as np
import os
import QD_functions as QD_scripts


def prep_core(core, core_folder, dummy=0, opt=True):
    """
    function that handles are core operations
    """    
    # Checks the if the dummy atom ligand placeholder is provided by its atomic number (int) or atomic symbol (string)
    # Returns an error if neither an integer nor string is provided
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum
    
    if opt:
        core = QD_scripts.global_minimum(core, core_folder)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core 
    # An additional dummy atom is added at the core center of mass for the purpose of orietating the ligands in QD_scripts.find_substructure(...)
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
            ligand = QD_scripts.global_minimum(ligand, ligand_folder)
        database_entry = QD_scripts.create_entry(ligand, opt)
        
    # find functional groups
    ligand.add_atom(Atom(atnum=0, coords=(ligand.get_center_of_mass())))
    ligand_list = QD_scripts.find_substructure(ligand)
    ligand_list = list(ligand_list)
    ligand_list.append(database_entry)
    
    return ligand_list


def prep_QD(core, ligand, core_indices, ligand_index, QD_folder):
    """
    function that handles all core+ligand operations
    """
    # Rotate and translate all ligands to their position on the core
    # Returns a list with sublist [0] containing the rotated ligands (PLAMS Molecules) and [1] the heteroatoms (PLAMS Atoms) of the rotated ligands to be attached to the core
    # All core dummy atoms are deleted
    core = copy.deepcopy(core)
    ligand_list = [QD_scripts.rotate_ligand(core, ligand, core_index, ligand_index) for core_index in core_indices]
    ligand_list = np.array(ligand_list).T.tolist()

    core.delete_atom(core[-1])

    # Prepare the .pdb filename (string)
    core_name = core.get_formula()
    ligand_name = ligand_list[0][0].get_formula() + '_@_' + ligand[ligand_index + 1].symbol + str(ligand_index + 1)
    pdb_name = str('core_' + core_name + '___ligand_' + ligand_name)

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule)
    QD = QD_scripts.combine_QD(core, ligand_list[0])

    # indices of all the atoms in the core and the ligand heteroatom anchor
    QD_indices = [QD.atoms.index(atom) + 1 for atom in ligand_list[1]]
    QD_indices += [i + 1 for i,atom in enumerate(core)]

    molkit.writepdb(QD, os.path.join(QD_folder, pdb_name + '.pdb'))
    print('core + ligands:\t\t\t' + pdb_name + '.pdb')   

    return QD, pdb_name, QD_indices





# The start
time_start = time.time()
print('\n')

# managing the result directories
dir_name_list = ['core', 'ligand', 'QD']
core_folder, ligand_folder, QD_folder = [QD_scripts.create_dir(name, path=r'/Users/basvanbeek/Documents/CdSe/Week_5') 
                                        for name in dir_name_list]

# Accepted inputs: .xyz/.pdb file, SMILES string, plain text file with SMILES strings or a list of aforementioned objects
input_cores = ['Cd68Se55.xyz']
input_ligands = 'OCCCCCCCCCO'

# Imports the cores and ligands
core_list = QD_scripts.read_mol(core_folder, input_cores)
ligand_list = QD_scripts.read_mol(ligand_folder, input_ligands, smiles_column=0, smiles_extension='.txt')

# Return the indices of the core dummy atoms
core_indices = [prep_core(core, core_folder, dummy='Cl', opt=False) for core in core_list]

# Open the ligand database and check if the specified ligand(s) is already present
database = QD_scripts.read_database(ligand_folder)
ligand_list = [prep_ligand(ligand, ligand_folder, database, opt=True) for ligand in ligand_list]

# formating of ligand_list
database_entries = [item[2] for item in ligand_list]
ligand_indices = [item[1] for item in ligand_list]
ligand_indices = list(itertools.chain(*ligand_indices))
ligand_list = [item[0] for item in ligand_list]
ligand_list = list(itertools.chain(*ligand_list))

# writing new entries to the ligand database
QD_scripts.write_database(database_entries, ligand_folder, database)

# combine the core with the ligands, yielding QD
QD_list = [prep_QD(core, ligand, core_indices[i], ligand_indices[j], QD_folder) for 
           i,core in enumerate(core_list) for j,ligand in enumerate(ligand_list)]

# formating of QD_list
QD_indices = [item[2] for item in QD_list]
pdb_name_list = [item[1] for item in QD_list]
QD_list = [item[0] for item in QD_list]

# optimize QD with the core frozen
#[QD_scripts.run_ams_job(QD, pdb_name_list[i], QD_folder, QD_indices[i]) for i,QD in enumerate(QD_list)]

# The End
time_end = time.time()
print("\nTotal elapsed time:\t\t" + "%.4f" % (time_end - time_start) + ' sec')
