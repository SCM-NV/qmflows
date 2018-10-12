import copy
import os
import itertools
import shutil
from functools import partial
import numpy as np

from scm.plams import (Molecule, Settings, AMSJob, init, finish)
from qmflows import molkit
from rdkit import Chem
from rdkit.Chem import Bond

import QD_database as qd_database

def remove_h(mol):
    """
    Remove all hydrogen atoms
    """
    for atom in mol:
        if atom.atnum == 1:
            mol.delete_atom(atom)

    return mol


def add_h(mol):
    """
    Add hydrogen atoms
    """
    mol = Chem.AddHs(molkit.to_rdmol(mol), addCoords=True)

    return molkit.from_rdmol(mol)


def concatenate_dict(dic):
    """
    Concatenates a list of dictionaries.
    """
    concact_dic = {}
    for item in dic:
        concact_dic.update(item)
    return concact_dic


def create_dir(dir_name, path=os.getcwd()):
    """
    Creates a new directory if this directory does not yet exist.
    """
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    return dir_path


def read_mol(input_mol, folder_path, column=0, row=0):
    """
    Checks by which means the input is provided and converts it into a molecule.
    Returns a dictionary of PLAMS molecules.
    """
    # Creates a dictionary of file extensions
    extension_dict = {'xyz' : read_mol_xyz, 'pdb' : read_mol_pdb, 'mol' : read_mol_mol,
                      'smiles' : read_mol_smiles, 'folder' : read_mol_folder,
                      'txt' : partial(read_mol_txt, row=row, column=column)}

    # Reads the input molecule(s), the method depending on the nature of the file extension
    # Returns a list of dictionaries
    mol_list = [extension_dict[input_mol[mol]](mol, folder_path) for mol in input_mol if
                input_mol[mol] in extension_dict]

    return concatenate_dict(mol_list)


def read_mol_xyz(mol, folder_path):
    """
    Read an .xyz file
    """
    folder_path = os.path.join(folder_path, mol)
    return {mol.rsplit('.', 1)[0]: Molecule(folder_path)}

def read_mol_pdb(mol, folder_path):
    """
    Read a .pdb file
    """
    folder_path = os.path.join(folder_path, mol)
    return {mol.rsplit('.', 1)[0]: molkit.readpdb(folder_path)}

def read_mol_mol(mol, folder_path):
    """
    Read a .mol file
    """
    folder_path = os.path.join(folder_path, mol)
    return {mol.rsplit('.', 1)[0]: molkit.from_rdmol(Chem.MolFromMolFile(folder_path))}

def read_mol_smiles(mol, folder_path):
    """
    Read a SMILES string
    """
    return {mol: molkit.from_smiles(mol)}

def read_mol_folder(mol, folder_path):
    """
    Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder
    """
    folder_path = os.path.join(folder_path, mol)
    file_dict = [{file.rsplit('.', 1)[0]: file.rsplit('.', 1)[1]} for
                  file in os.listdir(folder_path)]
    file_dict = concatenate_dict(file_dict)
    return read_mol(file_dict, folder_path)

def read_mol_txt(mol, folder_path, row, column):
    """
    Read a plain text file containing one or more SMILES strings
    """
    folder_path = os.path.join(folder_path, mol)
    with open(folder_path, 'r') as file:
        smiles_list = file.read().splitlines()
    smiles_list = [smiles.split()[column] for smiles in smiles_list[row:] if smiles]
    mol_list = [molkit.from_smiles(mol) for mol in smiles_list]
    smiles_list = [smile.replace('/', '') for smile in smiles_list]
    smiles_list = [smile.replace('\\', '') for smile in smiles_list]
    mol_dict = dict(zip(smiles_list, mol_list))
    return mol_dict


def set_pdb(mol, residue_name='RES', is_core=True):
    """
    Set a number of atomic properties
    """
    # Prepare a list of letters for pdb_info.Name
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet = [i + j for i in alphabet for j in alphabet]

    # Create a dictionary of elements and their formal atomic charge
    elements_dict = set_pdb_dict()

    # Set the atomic properties
    for i, atom in enumerate(mol):
        symbol = atom.symbol + alphabet[i] + '  '

        # Add a number of properties to atom
        atom.properties.pdb_info.ResidueName = residue_name
        atom.properties.pdb_info.Occupancy = 1.0
        atom.properties.pdb_info.TempFactor = 0.0
        atom.properties.pdb_info.ResidueNumber = 1
        atom.properties.pdb_info.Name = symbol[:4]
        atom.properties.pdb_info.ChainId = 'A'

        # Changes hydrogen and carbon from heteroatom to atom
        if atom.symbol == 'H' or atom.symbol == 'C':
            atom.properties.pdb_info.IsHeteroAtom = False

        # Sets the formal atomic charge
        if is_core and atom.atnum != 0:
            atom.properties.charge = elements_dict[atom.symbol]

    return mol


def set_pdb_dict():
    """
    Create a dictionary of elements and their formal atomic charge
    """
    # Create a list of atomic charges and elements
    charges = [1, 2, -3, -2, -1, 2]
    a = ['Li', 'Na', 'K', 'Rb', 'Cs']       # Alkaline metals
    b = ['Be', 'Mg', 'Ca', 'Sr', 'Ba']      # Alkaline earth metals
    c = ['N', 'P', 'As', 'Sb', 'Bi']        # Pnictogens
    d = ['O', 'S', 'Se', 'Te', 'Po']        # Chalcogens
    e = ['H', 'F', 'Cl', 'Br', 'I', 'At']   # Halogens
    f = ['Cd', 'Pb']                        # Misc

    # Combine the elements and atomic charges into a dictionary
    elements = [a, b, c, d, e, f]
    charges = [charges[i] for i, column in enumerate(elements) for element in column]
    elements = list(itertools.chain(*elements))

    return dict(zip(elements, charges))


def manage_ligand(ligand, ligand_name, ligand_folder, opt, database):
    """
    Pull the structure if a match has been found or alternatively optimize a new geometry.
    """
    # Searches for matches between the input ligand and the database; imports the structure

    ligand, match, pdb = qd_database.compare(ligand, ligand_name, ligand_folder, database)
    ligand_name = 'Ligand_' + ligand_name

    if ligand.properties.smiles:
        ligand_smiles = ligand.properties.smiles
    else:
        ligand_smiles = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(ligand)))

    # Optimize the ligand if no match has been found with the database
    if not match or not pdb:
        # Export the unoptimized ligand to a .pdb and .xyz file
        export_mol(ligand, ligand_folder, ligand_name, message='Ligand:\t\t\t\t')

        # If ligand optimization is enabled: Optimize the ligand, set pdb_info and export the result
        if opt:
            ligand = remove_h(ligand)
            ligand = molkit.global_minimum(ligand, n_scans=1, no_h=True)
            ligand = add_h(ligand)
            ligand = molkit.global_minimum(ligand, n_scans=1, no_h=True)
            set_pdb(ligand, residue_name='LIG', is_core=False)
            export_mol(ligand, ligand_folder, ligand_name + '.opt', message='Optimized ligand:\t\t')

        # Create an entry for in the database if no previous entries are present
        # or prints a warning if a structure is present in the database but the .pdb file is missing
        if not match and not pdb:
            entry = qd_database.entry(ligand, ligand_name, ligand_smiles, opt)
        else:
            entry = False
            print('\ndatabase entry exists for ' + ligand_name +
                  ' yet the corresponding .pdb file is absent. The geometry has been reoptimized.')
    else:
        entry = False

    return ligand, entry


def export_mol(mol, mol_folder, mol_name, message='Mol:\t\t\t\t'):
        """
        Write results to a .pdb and .xyz file
        """
        molkit.writepdb(mol, os.path.join(mol_folder, mol_name + '.pdb'))
        mol.write(os.path.join(mol_folder, mol_name + '.xyz'))
        print(str(message) + str(mol_name) + '.pdb')


def find_substructure(ligand, ligand_name, split=True):
    """
    Identify the ligand functional groups.
    """
    ligand_rdkit = molkit.to_rdmol(ligand)

    # Creates a list containing predefined functional groups, each saved as an rdkit molecule
    # IMPORTANT: The first atom should ALWAYS be the atom that should attach to the core
    if split:
        functional_group_list = ['[N+].[-]',
                                 'O[H]',
                                 'S[H]',
                                 'N[H]',
                                 'P[H]',
                                 '[O-].[+]',
                                 '[S-].[+]',
                                 '[N-].[+]',
                                 '[P-].[+]']
    else:
        functional_group_list = ['[N+]',
                                 'O[H]',
                                 'S[H]',
                                 'N[H]',
                                 'P[H]',
                                 '[O-]',
                                 '[S-]',
                                 '[N-]',
                                 '[P-]']

    functional_group_list = [Chem.MolFromSmarts(smarts) for smarts in functional_group_list]

    # Searches for functional groups (defined by functional_group_list) within the ligand
    # Duplicates are removed
    get_match = ligand_rdkit.GetSubstructMatches
    matches = [get_match(mol) for mol in functional_group_list]
    matches = list(itertools.chain(*matches))

    # Remove all duplicate matches
    ligand_indices = []
    for match in matches:
        if match[0] not in [item[0] for item in ligand_indices]:
            ligand_indices.append(match)

    if ligand_indices:
        ligand_list = [copy.deepcopy(ligand) for match in ligand_indices]
        ligand_list = [find_substructure_split(ligand, ligand_name, ligand_indices[i], split) for
                       i, ligand in enumerate(ligand_list)]
        ligand_list, ligand_names, ligand_indices = zip(*ligand_list)
        ligand_dict = dict(zip(ligand_names, ligand_list))
    else:
        print('No functional groups were found for ' + str(ligand.get_formula()))
        ligand_dict = {}
        ligand_indices = []

    return ligand_dict, ligand_indices


def find_substructure_split(ligand, ligand_name, ligand_index, split=True):
    """
    Delete the hydrogen or mono-/polyatomic counterion attached to the functional group
    Sets the charge of the remaining heteroatom to -1 if split=True
    """
    at1 = ligand[ligand_index[0] + 1]
    at2 = ligand[ligand_index[1] + 1]
    ligand_name = 'Ligand_' + ligand_name + '@' + at1.symbol + str(ligand_index[0] + 1)
    
    if split:
        if len(ligand.separate()) == 1:
            ligand.delete_atom(at2)
        else:
            mol1, mol2 = ligand.separate()
            if str(at1) in [str(atom) for atom in mol1]:
                ligand = mol1
            else:
                ligand = mol2
                
        # Check if the ligand heteroatom has a charge assigned, assigns a charge if not
        if not at1.properties.charge or at1.properties.charge == 0:
            at1.properties.charge = -1

    # Update the index of the ligand heteroatom
    ligand_atoms = [str(atom) for atom in ligand]
    ligand_index = ligand_atoms.index(str(at1)) + 1
    
    return ligand, ligand_name, ligand_index


def rotate_ligand(core, ligand, core_index, ligand_index, i):
    """
    Connects two molecules by alligning the vectors of two bonds.
    """
    ligand = copy.deepcopy(ligand)

    # Defines first atom on coordinate list (hydrogen),
    # The atom connected to it and vector representing bond between them
    core_at1 = core[core_index]         # core dummy atom
    core_at2 = core[-1]                 # core center of mass
    core_vector = core_at1.vector_to(core_at2)
    lig_at1 = ligand[ligand_index]  # ligand heteroatom
    lig_at2 = ligand[-1]                # ligand center of mass
    lig_vector = lig_at2.vector_to(lig_at1)

    # Rotation of ligand - aligning the ligand and core vectors
    rotmat = rotate_ligand_rotation(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(lig_at1.vector_to(core_at1))

    # Translation of the ligand
    hc_vec = lig_at1.vector_to(core_at1)
    ligand.translate(hc_vec)

    # Update the residue numbers
    for atom in ligand:
        atom.properties.pdb_info.ResidueNumber = i + 2

    # Deletes the core dummy atom and ligand center of mass
    ligand.delete_atom(lig_at2)
    core.delete_atom(core_at1)

    return ligand, lig_at1


def rotate_ligand_rotation(vec1, vec2):
    """
    Calculates the rotation matrix rotating *vec1* to *vec2*.
    Vectors can be any containers with 3 numerical values. They don't need to be normalized.
    Returns 3x3 numpy array.
    """
    a = np.array(vec1) / np.linalg.norm(vec1)
    b = np.array(vec2) / np.linalg.norm(vec2)
    v1, v2, v3 = np.cross(a, b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])

    return np.identity(3) + M + np.dot(M, M)/(1+np.dot(a, b))


def combine_qd(core, ligand_list):
    """
    Combine the rotated ligands with the core, creating a bond bewteen the core and ligand.
    """
    qd = copy.deepcopy(core)

    # Create a list of ligand atoms and intraligand bonds
    ligand_bonds = np.concatenate([ligand.bonds for ligand in ligand_list])
    ligand_atoms = np.concatenate(ligand_list)

    # Combined the ligand bond and atom list with the core
    for atom in ligand_atoms:
        qd.add_atom(atom)
    for bond in ligand_bonds:
        qd.add_bond(bond)

    return qd


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
        print('One or more ADF environment variables have not been set, aborting ' +
              'geometry optimization.')
        return False
    else:
        return True


def ams_job(qd, qd_name, qd_folder, qd_indices, maxiter=1000):
    """
    Converts PLAMS connectivity into adf .run script connectivity.
    """
    # Create a list of aromatic bond indices
    qd_rdkit = molkit.to_rdmol(qd)
    aromatic = [Bond.GetIsAromatic(bond) for bond in qd_rdkit.GetBonds()]
    aromatic = [i for i, item in enumerate(aromatic) if item]

    # Create a connectivity list; aromatic bonds get a bond order of 1.5
    at1 = [qd.atoms.index(bond.atom1) + 1 for bond in qd.bonds]
    at2 = [qd.atoms.index(bond.atom2) + 1 for bond in qd.bonds]
    bonds = [bond.order for bond in qd.bonds]
    for i, bond in enumerate(qd.bonds):
        if i in aromatic:
            bonds[i] = 1.5
    bonds = [str(at1[i]) + ' ' + str(at2[i]) + ' ' + str(bond) for i, bond in enumerate(bonds)]

    # Launch the AMS UFF constrained geometry optimization
    output_mol = ams_job_run(qd, qd_name, qd_folder, qd_indices, bonds, maxiter)

    # Update the atomic coordinates of qd
    for i, atom in enumerate(qd):
        atom.move_to(output_mol[i + 1])

    # Write the reuslts to an .xyz and .pdb file
    export_mol(qd, qd_folder, qd_name, message='Optimized core + ligands:\t\t')

    return qd


def ams_job_run(qd, pdb_name, qd_folder, qd_indices, bonds, maxiter):
    """
    Runs the AMS UFF constrained geometry optimization.
    """
    # General AMS settings
    s = Settings()
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.Constraints.Atom = qd_indices
    s.input.ams.System.BondOrders._1 = bonds
    s.input.ams.GeometryOptimization.MaxIterations = maxiter
    s.input.ams.Properties.Gradients = 'Yes'

    # Settings specific to UFF
    s.input.uff.Library = 'UFF'

    # Run the job
    init(path=qd_folder, folder=pdb_name)
    job = AMSJob(molecule=qd, settings=s, name=pdb_name)
    results = job.run()
    output_mol = results.get_main_molecule()
    finish()

    # Copy the resulting .rkf and .out files and delete the PLAMS directory
    shutil.copy2(results['ams.rkf'], os.path.join(qd_folder, pdb_name + '.opt.rkf'))
    shutil.copy2(results[pdb_name + '.out'], os.path.join(qd_folder, pdb_name + '.opt.out'))
    shutil.rmtree(os.path.join(qd_folder, pdb_name))

    return output_mol
