import copy
import os
import itertools
import time
import string

from scm.plams import *
from qmflows import molkit

import QD_functions as qd_scripts


def prep_core(core, core_folder, dummy=0, opt=True):
    """
    Function that handles all core operations.
    """
    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum

    # Optimize the core with RDKit UFF if opt = True. Returns a RDKit molecule
    if opt:
        core = qd_scripts.global_minimum(core, core_folder)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    core_indices = [(i + 1) for i, atom in enumerate(core.atoms) if atom.atnum == dummy]
    core_indices.reverse()
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Sets core pdb_info 
    for atom in core:
        atom.properties.pdb_info.ResidueName = 'COR'
        atom.properties.pdb_info.ResidueNumber = 1
        atom.properties.pdb_info.Occupancy = 1.0
        atom.properties.pdb_info.ChainId = 'A'
        name = atom.symbol + '   '
        atom.properties.pdb_info.Name = name[:4]

    # Sets formal atomic charges
    for atom in core:
        if atom.symbol in ['Li', 'Na', 'K', 'Rb', 'Cs']:
            atom.properties.charge = 1
        elif atom.symbol in ['Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Cd', 'Pb']:
            atom.properties.charge = 2
        elif atom.symbol in ['N', 'P', 'As', 'Sb', 'Bi']:
            atom.properties.charge = -3
        elif atom.symbol in ['O', 'S', 'Se', 'Te', 'Po']:
            atom.properties.charge = -2
        elif atom.symbol in ['H', 'F', 'Cl', 'Br', 'I', 'At']:
            atom.properties.charge = -1

    # Returns an error if no dummy atoms were found
    if not core_indices:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')
    else:
        return core_indices


def prep_ligand(ligand, ligand_folder, database, opt=True):
    """
    Function that handles all ligand operations,
    """
    # Checks if the database exists and if the ligand is already present.
    # Returns a list of booleans.
    if database:
        matches = [molkit.to_rdmol(ligand).HasSubstructMatch(mol) for mol in database[4]]
    else:
        matches = [False]

    # If the ligand is already present in the database: append the database.
    # If not, optimize the ligand and create a new entry for the database.
    if any(matches):
        index = matches.index(True)
        ligand = molkit.readpdb(os.path.join(ligand_folder, str(database[3][index])))
        database_entry = False
    else:
        # Export the unoptimized ligand to a .pdb and .xyz file
        ligand_name = 'ligand_' + ligand.get_formula()
        molkit.writepdb(ligand, os.path.join(ligand_folder, ligand_name + '.pdb'))
        print('Ligand:\t\t\t\t' + str(ligand_name) + '.pdb')
        if opt:
            # Optimize the ligand
            ligand = qd_scripts.global_minimum(ligand, ligand_folder)

            # Set pdb_info
            alphabet = list(string.ascii_lowercase)
            alphabet = [i + j for i in alphabet for j in alphabet]
            for i, atom in enumerate(ligand):
                atom.properties.pdb_info.ResidueName = 'LIG'
                atom.properties.pdb_info.Occupancy = 1.0
                atom.properties.pdb_info.ResidueNumber = 1
                symbol = atom.symbol + alphabet[i] + '  '
                atom.properties.pdb_info.Name = symbol[:4]
                if atom.symbol == 'H' or atom.symbol == 'C':
                    atom.properties.pdb_info.IsHeteroAtom = False
                atom.properties.pdb_info.ChainId = 'A'
                
            # Export the optimized ligand to a .pdb and .xyz file
            molkit.writepdb(ligand, os.path.join(ligand_folder, ligand_name + '.opt.pdb'))
            print('Optimized ligand:\t\t' + str(ligand_name) + '.opt.pdb')
        database_entry = qd_scripts.create_entry(ligand, opt)



    # Identify functional groups within the ligand.
    ligand.add_atom(Atom(atnum=0, coords=ligand.get_center_of_mass()))
    ligand_list = qd_scripts.find_substructure(ligand)
    ligand_list = list(ligand_list)
    ligand_list.append(database_entry)


    return ligand_list


def prep_qd(core, ligand, core_indices, ligand_index, qd_folder):
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
    core_name = core.get_formula()
    ligand_formula = ligand_list[0].get_formula()
    ligand_heteroatom = ligand[ligand_index + 1].symbol
    ligand_name = ligand_formula + '_@_' + ligand_heteroatom + str(ligand_index + 1)
    pdb_name = str('core_' + core_name + '___ligand_' + ligand_name)

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    qd = qd_scripts.combine_qd(core, ligand_list)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in ligand_indices]
    qd_indices += [i + 1 for i, atom in enumerate(core)]

    molkit.writepdb(qd, os.path.join(qd_folder, pdb_name + '.pdb'))
    print('core + ligands:\t\t\t' + pdb_name + '.pdb')

    return qd, pdb_name, qd_indices


def prep_prep(path, dir_name_list, input_cores, input_ligands, smiles_extension, smiles_column,
              dummy, core_opt, ligand_opt, qd_opt, maxiter):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.
    """
    # The start
    time_start = time.time()
    print('\n')

    # Managing the result directories
    core_folder, ligand_folder, qd_folder = [qd_scripts.create_dir(name, path=path)
                                             for name in dir_name_list]

    # Imports the cores and ligands
    core_list = qd_scripts.read_mol(core_folder, input_cores)
    ligand_list = qd_scripts.read_mol(ligand_folder, input_ligands, smiles_column, smiles_extension)

    # Return the indices of the core dummy atoms
    core_indices = [prep_core(core, core_folder, dummy, core_opt) for core in core_list]

    # Open the ligand database and check if the specified ligand(s) is already present
    database = qd_scripts.read_database(ligand_folder)
    ligand_list = [prep_ligand(ligand, ligand_folder, database, ligand_opt) for
                   ligand in ligand_list]

    # Formating of ligand_list
    ligand_list, ligand_indices, database_entries = zip(*ligand_list)
    ligand_indices = list(itertools.chain(*ligand_indices))
    ligand_list = itertools.chain(*ligand_list)

    # Write new entries to the ligand database
    qd_scripts.write_database(database_entries, ligand_folder, database)

    # Combine the core with the ligands, yielding qd
    qd_list = [prep_qd(core, ligand, core_indices[i], ligand_indices[j], qd_folder) for
               i, core in enumerate(core_list) for j, ligand in enumerate(ligand_list)]

    # Formating of qd_list
    qd_list, pdb_name_list, qd_indices = zip(*qd_list)

    # Optimize qd with the core frozen
    for i, qd in enumerate(qd_list):
        qd_scripts.run_ams_job(qd, pdb_name_list[i], qd_folder, qd_indices[i], maxiter, qd_opt)

    # The End
    time_end = time.time()
    print('\nTotal elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')




"""
path =              The path where the input and output directories will be saved.
dir_name_list =     Names of the to be created directories in path. Set to os.getcwd() to use the
                    current directory.
input_cores =       The input core(s) as either .xyz, .pdb, .mol, SMILES string, plain text file 
                    with SMILES string or a list containing any of the above objects.
input_ligands =     Same as input_cores, except for the ligand(s).
smiles_extension =  Extension of a SMILES string containg plain text file. Relevant if such a file 
                    is chosen for input_cores or input_ligands.
smiles_column =     The column containing the SMILES string in a plain text file. Relevant if such 
                    a file is chosen for input_cores or input_ligands.
dummy =             The atomic number of atomic symbol of the atoms in the core that should be
                    should be replaced with ligands.
core_opt =          Attempt to find the core global minimum using RDKit UFF. 
                    WARNING: enabling this will probably ruin the core if care is not taken.
                    Should work fine for organic cores.
ligand_opt =        Attempt to find the ligand global minimum using RDKit UFF.
qd_opt =            Optimize the quantum dot (qd, i.e core + all ligands) using ADF UFF.
maxiter =           The maximum number of geometry iteration for qd_opt.
"""

# Argument list
path = r'/Users/basvanbeek/Documents/CdSe/Week_5'
dir_name_list = ['core', 'ligand', 'QD']
input_cores = 'Cd68Se55.xyz'
input_ligands = 'OCCCCCCCCC'
smiles_extension = '.txt'
smiles_column = 0
dummy = 'Cl'
core_opt = False
ligand_opt = True
qd_opt = False
maxiter = 100

# Runs the script: add ligand to core and optimize (UFF) the resulting qd with the core frozen
prep_prep(path, dir_name_list, input_cores, input_ligands, smiles_extension, smiles_column, dummy,
          core_opt, ligand_opt, qd_opt, maxiter)
