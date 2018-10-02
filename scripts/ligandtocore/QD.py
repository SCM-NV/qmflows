from scm.plams import *
from qmflows import molkit
import time
import copy
import numpy as np
import os
import QD_functions as QD


def prep_core(core, ligand_list, ligand_folder, core_folder, core_ligand_folder, dummy=0, ligand_opt=True, core_opt=False, core_ligand_opt=True):
    """
    Identify the ligand anchoring sites on the core, optimize the ligands, find the anchoring sites on the ligand and finnally attach the ligands to the core
    """    
    # Checks the if the dummy atom ligand placeholder is provided by its atomic number (int) or atomic symbol (string)
    # Returns an error if neither an integer nor string is provided
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum
    
    # Returns the indices (integer) of all dummy atom ligand placeholders in the core 
    # An additional dummy atom is added at the core center of mass for the purpose of orietating the ligands in QD.find_substructure(...)
    core_indices = [(i + 1) for i, atom in enumerate(core.atoms) if atom.atnum == dummy]
    core_indices.reverse()
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # checks if the ligand is already present in the database, either creating a new entry or appending an existing entry
    ligand_list = copy.deepcopy(ligand_list)
    ligand_list = QD.check_database(ligand_list, ligand_folder, ligand_opt, database_name='ligand_database.txt')

    # Searches for the global minimum (using UFF) by systematically evaluating all dihedral angles
    # Returns the optimized structure (PLAMS Molecule)
    # Previously optimized structures are not reoptimized
    if ligand_opt:
        opt = [QD.global_minimum(ligand[0], ligand_folder) for ligand in ligand_list if ligand[1]]
        no_opt = [ligand[0] for ligand in ligand_list if not ligand[1]] 
        ligand_list = opt + no_opt
    if core_opt:
        core = QD.global_minimum(core, core_folder)
    
    # Identify the functional groups within the ligand that can bond with the core
    # One additional copy of the ligand is added to ligand_list for each functional group beyond the first
    # returns a list with sublist [0] containing the ligand (PLAMS Molecule) and [1, 2 & 3] the indices (integer) of the functional group hydrogen, heteroatom and beta carbon
    # Returns an error if no functional groups are found for a specific ligand
    ligand_list = [QD.find_substructure(ligand) for ligand in ligand_list if QD.find_substructure(ligand)]
    ligand_list = np.concatenate(ligand_list)

    # Rotate the ligands and combine them with the core, returing the resulting structure (PLAMS Molecule)
    core_ligand_list = [prep_ligand(core, ligand, core_indices, core_ligand_folder, core_ligand_opt) for ligand in ligand_list]

    return core_ligand_list


def prep_ligand(core, ligand, core_indices, core_ligand_folder, core_ligand_opt):
    """
    add all ligands to the core
    """
    # Rotate and translate all ligands to their position on the core
    # Returns a list with sublist [0] containing the rotated ligands (PLAMS Molecules) and [1] the heteroatoms (PLAMS Atoms) of the rotated ligands to be attached to the core
    # All core dummy atoms are deleted
    core = copy.deepcopy(core)
    ligand[0].add_atom(Atom(atnum=0, coords=(ligand[0].get_center_of_mass())))
    ligand_list = [QD.rotate_ligand(core, ligand, index) for index in core_indices]
    ligand_list = np.array(ligand_list).T.tolist()
    core.delete_atom(core[-1])

    # Prepare the .pdb filename (string)
    core_name = core.get_formula()
    ligand_name = ligand_list[0][0].get_formula() + '_@_' + ligand[0][ligand[2]].symbol + str(ligand[2])
    pdb_name = str('core_' + core_name + '__&__ligand_' + ligand_name)

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule)
    core_ligand = QD.combine_core_ligand(core, ligand_list[0])
    
    # Create list with sublists containg [0] ligand heteroatom indices (integer), [1] core selenium indices (integer) and [2] core cadmium indices (integer)
    atnum_list = []
    for atom in core:
        if atom.atnum not in atnum_list:
            atnum_list.append(atom.atnum)

    core_ligand_indices = [[core_ligand.atoms.index(atom) for atom in ligand_list[1]]]
    core_ligand_indices += [[i for i,atom in enumerate(core) if atom.atnum == atnum] for atnum in atnum_list]
    
    # Assign residue names and formal atomic charges
    #core_ligand = molkit.to_rdmol(core_ligand)
    #core_ligand = QD.prepare_pdb(core_ligand, core, ligand, core_ligand_indices)
    
    # Optimize the combined core and ligands, writing the results to a .pdb file
    #molkit.writepdb(core_ligand, os.path.join(core_ligand_folder, pdb_name + '.pdb'))
    print('core + ligands:\t\t\t' + pdb_name + '.pdb')
    #core_ligand = molkit.from_rdmol(core_ligand)
    if core_ligand_opt:
        core_ligand = QD.optimize_core_ligand(core_ligand, core_ligand_indices, maxiter=200)
        molkit.writepdb(core_ligand, os.path.join(core_ligand_folder, pdb_name + '.opt.pdb'))
        print('\nOptimized core + ligands:\t' + pdb_name + '.opt.pdb')

    return core_ligand





# The start
time_start = time.time()
print('\n')

dir_name_list = ['ligand', 'core', 'core_ligand']
dir_path_list = [QD.create_dir(name, path=r'C:\Users\hardd\Documents\Computational\CdSe') for name in dir_name_list]

# Accepted inputs: .xyz/.pdb file, SMILES string, plain text file with SMILES strings or a list of aforementioned objects
input_ligands = 'OC(C1=CC=CC=C1)=O'
input_cores = ['cube.xyz', 'cube_Mg.xyz']

# Imports the cores and ligands
ligand_list = QD.read_mol(dir_path_list[0], input_ligands, smiles_column=0, smiles_extension='.txt')
core_list = QD.read_mol(dir_path_list[1], input_cores)

# Copies of ligands are added to copies of core molecules; this process is repeated until all dummy atoms are substituted for ligands
core_ligand_list = [prep_core(core, ligand_list, dir_path_list[0], dir_path_list[1], dir_path_list[2], dummy='Rb', ligand_opt=True, core_opt=False, core_ligand_opt=False) for core in core_list]
core_ligand_list = np.concatenate(core_ligand_list)
[QD.run_ams_job(core_ligand) for core_ligand in core_ligand_list]

# The End
time_end = time.time()
print("\nTotal elapsed time:\t\t" + "%.4f" % (time_end - time_start) + ' sec')
