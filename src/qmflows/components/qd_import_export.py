__all__ = ['read_mol', 'set_prop', 'create_dir']

import os
import itertools
import pandas as pd

from scm.plams import Molecule
from scm.plams.core import errors
import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem


def read_mol(input_mol, folder_path, is_core=False):
    """
    Checks the filetypes of the input molecules, sets their properties and
    returns a list of plams molecules.
    """
    # Creates a dictionary of file extensions
    extension_dict = {'xyz': read_mol_xyz,
                      'pdb': read_mol_pdb,
                      'mol': read_mol_mol,
                      'smiles': read_mol_smiles,
                      'folder': read_mol_folder,
                      'txt': read_mol_txt,
                      'plams_mol': read_mol_plams,
                      'rdmol': read_mol_rdkit}

    # Reads the input molecule(s), the method depending on the nature of the file extension
    # Returns a list of dictionaries
    input_mol = [read_mol_defaults(mol, folder_path, is_core) for mol in input_mol]
    input_mol = [read_mol_extension(mol) for mol in input_mol]
    input_mol = dict_concatenate(input_mol)

    # Create a list of PLAMS molecules
    mol_list = []
    for i, mol in enumerate(input_mol):
        mol_dict = input_mol[mol]
        # Check if the provided key is available in extension_dict
        try:
            read = extension_dict[mol_dict['file_type']]
        except KeyError as ex:
            print(str(type(ex).__name__) + ':\t' + str(ex) + '\n')
        # Convert mol into either a list or PLAMS molecule
        if read:
            mol = read(mol, mol_dict)
            # if mol is a PLAMS molecule
            if isinstance(mol, Molecule):
                # Guess bonds if guess_bonds=True
                if mol_dict['guess_bonds']:
                    mol.guess_bonds()
                # Guess bonbs if an xyz file was provided and the user did not specifiy guess_bonds
                if mol_dict['file_type'] == 'xyz':
                    if isinstance(mol_dict['user_args'], str):
                        mol.guess_bonds()
                    else:
                        keys = list(mol_dict['user_args'][1].keys())
                        if 'guess_bonds' not in keys:
                            mol.guess_bonds()
                set_prop(mol, mol_dict)
                mol_list.append(mol)
            # if mol is a list
            if isinstance(mol, list):
                mol_list += mol

    # Raises an error if mol_list is empty
    if not mol_list:
        core_ligand = {True: 'cores', False: 'ligands'}
        error = 'No valid input ' + core_ligand[is_core] + ' were found, aborting run'
        raise IndexError(error)

    return mol_list


def read_mol_defaults(mol, folder_path, is_core=False):
    """
    Assign all (default) arguments to mol
    """
    # A dictionary of default arguments
    kwarg = {'file_type': '',
             'mol_name': mol,
             'row': 0,
             'column': 0,
             'sheet_name': 'Sheet1',
             'guess_bonds': False,
             'is_core': is_core,
             'folder_path': folder_path,
             'mol_path': '',
             'user_args': mol}

    # Update the dictionary of default arguments based on used-specified arguments
    if isinstance(mol, (list, tuple)):
        kwarg.update(mol[1])
        mol = mol[0]
    if mol is None:
        mol = ''

    kwarg['mol_path'] = os.path.join(kwarg['folder_path'], str(mol))

    return {mol: kwarg}


def read_mol_extension(mol_dict):
    """
    Identifies the filetype of mol.
    """
    # Identify the filetype of mol_name
    mol = list(mol_dict.keys())[0]
    if os.path.isfile(mol_dict[mol]['mol_path']):
        file_type = mol.rsplit('.', 1)[-1]
        mol_name = mol.rsplit('.', 1)[0]
    elif os.path.isdir(mol_dict[mol]['mol_path']):
        file_type = 'folder'
        mol_name = mol
    elif isinstance(mol, Molecule):
        file_type = 'plams_mol'
        mol_name = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))
    elif isinstance(mol, Chem.rdchem.Mol):
        file_type = 'rdmol'
        mol_name = Chem.MolToSmiles(Chem.RemoveHs(mol))
    else:
        file_type = 'smiles'
        mol_name = smiles_name(mol)

    mol_dict[mol]['file_type'] = file_type
    mol_dict[mol]['mol_name'] = mol_name

    return mol_dict


def smiles_name(mol):
    """
    Turn a SMILES string into an acceptable filename.
    """
    mol_name = mol.replace('(', '[').replace(')', ']')
    cis_trans = [item for item in mol if item is '/' or item is '\\']
    if cis_trans:
        cis_trans = [item + cis_trans[i*2+1] for i, item in enumerate(cis_trans[::2])]
        cis_trans_dict = {'//': 'trans-', '/\\': 'cis-'}
        for item in cis_trans[::-1]:
            mol_name = cis_trans_dict[item] + mol_name
        mol_name = mol_name.replace('/', '').replace('\\', '')

    return mol_name


def read_mol_xyz(mol, mol_dict):
    """
    Read an .xyz file
    """
    try:
        return Molecule(mol_dict['mol_path'], inputformat='xyz')
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_xyz.__code__, ex, mol_dict, mol_dict['mol_path'])


def read_mol_pdb(mol, mol_dict):
    """
    Read a .pdb file
    """
    try:
        return molkit.readpdb(mol_dict['mol_path'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_pdb.__code__, ex, mol_dict, mol_dict['mol_path'])


def read_mol_mol(mol, mol_dict):
    """
    Read a .mol file
    """
    try:
        return molkit.from_rdmol(Chem.MolFromMolFile(mol_dict['mol_path'], removeHs=False))
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_mol.__code__, ex, mol_dict, mol_dict['mol_path'])


def read_mol_smiles(mol, mol_dict):
    """
    Read a SMILES string
    """
    try:
        return molkit.from_smiles(mol)
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_smiles.__code__, ex, mol_dict, mol)


def read_mol_plams(mol, mol_dict):
    """
    Read a PLAMS molecule
    """
    try:
        return mol
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_plams.__code__, ex, mol_dict, mol)


def read_mol_rdkit(mol, mol_dict):
    """
    Read a RDKit molecule
    """
    try:
        return molkit.from_rdmol(mol)
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_rdkit.__code__, ex, mol_dict, mol)


def read_mol_folder(mol, mol_dict):
    """
    Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder
    """
    try:
        file_list = [[file, mol_dict] for file in os.listdir(mol_dict['mol_path'])]
        return read_mol(file_list, mol_dict['folder_path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_folder.__code__, ex, mol_dict, mol_dict['mol_path'])


def read_mol_txt(mol, mol_dict):
    """
    Read a plain text file containing one or more SMILES strings
    """
    try:
        with open(mol_dict['mol_path'], 'r') as file:
            file_list = file.read().splitlines()
        file_list = [file.split()[mol_dict['column']] for file in file_list[mol_dict['row']:] if
                     file]
        file_list = [[file, mol_dict] for file in file_list if len(file) >= 2]
        return read_mol(file_list, mol_dict['folder_path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_txt.__code__, ex, mol_dict, mol_dict['mol_path'])


def read_mol_excel(mol, mol_dict):
    """
    Read a plain text file containing one or more SMILES strings
    """
    try:
        file_list = pd.read_excel(mol_dict['mol_path'], sheet_name=mol_dict['sheet_name'])
        file_list = [[file, mol_dict] for file in file_list[mol_dict['column']][mol_dict['row']:]]
        return read_mol(file_list, mol_dict['folder_path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_excel.__code__, ex, mol_dict, mol_dict['mol_path'])


def set_prop(mol, mol_dict):
    """
    Set molecular and atomic properties
    """
    if mol_dict.get('is_core'):
        residue_name = 'COR'
        mol.properties.name = 'Core_' + mol_dict['mol_name']
    else:
        residue_name = 'LIG'
        mol.properties.name = 'Ligand_' + mol_dict['mol_name']

    mol.properties.formula = mol.get_formula()
    mol.properties.source_folder = mol_dict['folder_path']
    if not mol.properties.smiles:
        mol.properties.smiles = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))

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


def print_exception(func, ex, mol_dict, name):
    """
    Manages the printing of exceptions upon failing to import a molecule.
    """
    extension_dict = {'read_mol_xyz': '.xyz file', 'read_mol_pdb': '.pdb file',
                      'read_mol_mol': '.mol file', 'read_mol_smiles': 'SMILES string',
                      'read_mol_folder': 'folder', 'read_mol_txt': '.txt file',
                      'read_mol_excel': '.xlsx file', 'read_mol_plams': 'PLAMS molecule',
                      'read_mol_rdkit': 'RDKit molecule'}
    print(str(type(ex).__name__), str(ex))
    print('function:', str(func.co_name) + str(func.co_varnames[:func.co_argcount]))
    # print('\t\targ1(mol):', type(mol), mol)
    # print('\t\targ2(mol_dict):', type(mol_dict), mol_dict)
    print('Warning:',  name, 'not recognized as a valid', extension_dict[func.co_name], '\n')
    return []


def dict_concatenate(dic):
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
    if not os.path.exists(path):
        error = path + ' not found, aborting run'
        raise FileNotFoundError(error)
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path
