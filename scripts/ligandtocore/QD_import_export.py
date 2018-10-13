import os
import itertools

from functools import partial
from scm.plams import Molecule
from qmflows import molkit
from rdkit import Chem


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


def read_mol(input_mol, folder_path, column=0, row=0, is_core=False):
    """
    Checks by which means the input is provided and converts it into a molecule.
    Returns a list of PLAMS molecules.
    """
    # Creates a dictionary of file extensions
    extension_dict = {'xyz': read_mol_xyz, 'pdb': read_mol_pdb, 'mol': read_mol_mol,
                      'smiles': read_mol_smiles, 'folder': read_mol_folder,
                      'txt': partial(read_mol_txt, row=row, column=column)}

    # Reads the input molecule(s), the method depending on the nature of the file extension
    # Returns a list of dictionaries
    mol_list = []
    for mol in input_mol:
        if isinstance(input_mol[mol], str):
            if input_mol[mol] in extension_dict:
                read = extension_dict[input_mol[mol]]
                mol_list.append(read(mol, folder_path, is_core))
        elif isinstance(input_mol[mol], list):
            if input_mol[mol][0] in extension_dict:
                read = extension_dict[input_mol[mol][0]]
                guess_bonds = input_mol[mol][1]
                mol_list.append(read(mol, folder_path, is_core, guess_bonds))
    return list(itertools.chain(*mol_list))


def read_mol_xyz(mol_name, folder_path, is_core=False, guess_bonds=True):
    """
    Read an .xyz file
    """
    mol_path = os.path.join(folder_path, mol_name)
    mol = Molecule(mol_path)
    mol_name = mol_name.rsplit('.', 1)[0]
    if guess_bonds:
        mol.guess_bonds()
    set_prop(mol, mol_name, folder_path, is_core)
    return [mol]


def read_mol_pdb(mol_name, folder_path, is_core=False, guess_bonds=False):
    """
    Read a .pdb file
    """
    mol_path = os.path.join(folder_path, mol_name)
    mol = molkit.readpdb(mol_path)
    mol_name = mol_name.rsplit('.', 1)[0]
    if guess_bonds:
        mol.guess_bonds()
    set_prop(mol, mol_name, folder_path, is_core)
    return [mol]


def read_mol_mol(mol_name, folder_path, is_core=False, guess_bonds=False):
    """
    Read a .mol file
    """
    mol_path = os.path.join(folder_path, mol_name)
    mol = molkit.from_rdmol(Chem.MolFromMolFile(mol_path))
    mol_name = mol_name.rsplit('.', 1)[0]
    if guess_bonds:
        mol.guess_bonds()
    set_prop(mol, mol_name, folder_path, is_core)
    return [mol]


def read_mol_smiles(mol_name, folder_path, is_core=False, guess_bonds=False):
    """
    Read a SMILES string
    """
    mol = molkit.from_smiles(mol_name)
    mol_name = mol_name.replace('/', '').replace('\\', '')
    if guess_bonds:
        mol.guess_bonds()
    set_prop(mol, mol_name, folder_path, is_core)
    return [mol]


def read_mol_folder(mol, folder_path, is_core=False):
    """
    Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder
    Will search in folder_path/mol/ if it exists or alternatively just in folder_path/
    """
    if os.path.isfile(os.path.join(folder_path, mol)):
        folder_path = os.path.join(folder_path, mol)
    file_dict = [{file: file.rsplit('.', 1)[1]} for file in os.listdir(folder_path)]
    file_dict = concatenate_dict(file_dict)
    return read_mol(file_dict, folder_path, is_core)


def read_mol_txt(mol, folder_path, is_core=False, row=0, column=0):
    """
    Read a plain text file containing one or more SMILES strings
    """
    mol_path = os.path.join(folder_path, mol)
    with open(mol_path, 'r') as file:
        smiles_list = file.read().splitlines()
    smiles_list = [smiles.split()[column] for smiles in smiles_list[row:] if smiles]
    mol_dict = [{mol: 'smiles'} for mol in smiles_list]
    mol_dict = concatenate_dict(mol_dict)
    return read_mol(mol_dict, folder_path, is_core)


def set_prop(mol, mol_name, folder_path, is_core=False):
    """
    Set a molecular and atomic properties
    """
    if is_core:
        residue_name = 'COR'
        mol.properties.name = 'Core_' + mol_name
    else:
        residue_name = 'LIG'
        mol.properties.name = 'Ligand_' + mol_name

    mol.properties.formula = mol.get_formula()
    mol.properties.source_folder = folder_path
    if not mol.properties.smiles:
        Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))

    # Prepare a list of letters for pdb_info.Name
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet = [i + j for i in alphabet for j in alphabet]

    # Create a dictionary of elements and their formal atomic charge
    elements_dict = set_prop_dict()

    # Set the atomic properties
    for i, atom in enumerate(mol):
        set_prop_atom(atom, alphabet[i], residue_name, elements_dict)


def set_prop_atom(atom, alphabet, residue_name, elements_dict):
    """
    Set atomic properties.
    """
    symbol = atom.symbol + alphabet + '  '

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
    if not atom.properties.charge:
        if atom.symbol in elements_dict:
            total_bonds = [bond.order for bond in atom.bonds]
            atom.properties.charge = elements_dict[atom.symbol] + int(sum(total_bonds))
        else:
            atom.properties.charge = 0


def set_prop_dict():
    """
    Create a dictionary of elements and their formal atomic charge.
    """
    # Create a list of atomic charges and elements
    charges = [1, 2, -3, -2, -1, 2]
    group01 = ['Li', 'Na', 'K', 'Rb', 'Cs']       # Alkaline metals
    group02 = ['Be', 'Mg', 'Ca', 'Sr', 'Ba']      # Alkaline earth metals
    group15 = ['N', 'P', 'As', 'Sb', 'Bi']        # Pnictogens
    group16 = ['O', 'S', 'Se', 'Te', 'Po']        # Chalcogens
    group17 = ['H', 'F', 'Cl', 'Br', 'I', 'At']   # Halogens
    misc = ['Cd', 'Pb']                        # Misc

    # Combine the elements and atomic charges into a dictionary
    elements = [group01, group02, group15, group16, group17, misc]
    charges = [charges[i] for i, column in enumerate(elements) for element in column]
    elements = list(itertools.chain(*elements))

    return dict(zip(elements, charges))


def export_mol(mol, message='Mol:\t\t\t\t'):
    """
    Write results to a .pdb and .xyz file
    """
    mol_name = mol.properties.name
    mol_path = os.path.join(mol.properties.source_folder, mol_name)
    molkit.writepdb(mol, mol_path + '.pdb')
    mol.write(mol_path + '.xyz')
    print(str(message) + str(mol_name) + '.pdb')
