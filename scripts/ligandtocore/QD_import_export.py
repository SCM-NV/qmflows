import os
import itertools
import pandas as pd

from scm.plams import Molecule
from scm.plams.core import errors
from qmflows import molkit
from rdkit import Chem


def print_exception(func, ex, mol_path, kwargs):
    extension_dict = {'read_mol_xyz': '.xyz file', 'read_mol_pdb': '.pdb file',
                      'read_mol_mol': '.mol file', 'read_mol_smiles': 'SMILES string',
                      'read_mol_folder': 'folder', 'read_mol_txt': '.txt file',
                      'read_mol_excel': '.xlsx file', 'read_mol_plams': 'PLAMS molecule',
                      'read_mol_rdkit': 'RDKit molecule'}
    print(str(type(ex).__name__), str(ex))
    print('\tfunction:', str(func.co_name) + str(func.co_varnames[:func.co_argcount]))
    for i, item in enumerate(func.co_varnames[:func.co_argcount]):
        print('\t\targ' + str(i) + '(' + func.co_varnames[i] + ')' + ':', type(kwargs[i]), kwargs[i])
    print('Warning:', mol_path, 'is not recognize as a valid', extension_dict[func.co_name], '\n')
    return []


def dict_concatenate(dic):
    """
    Concatenates a list of dictionaries.
    """
    concact_dic = {}
    for item in dic:
        concact_dic.update(item)
    return concact_dic


def dict_defaults(dic, dic_ref):
    for item in dic_ref:
        if not dic.get(item):
            dic.update({item: dic_ref[item]})
    return dic


def create_dir(dir_name, path=os.getcwd()):
    """
    Creates a new directory if this directory does not yet exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path, 'not found, aborting run')
    dir_path = os.path.join(path, str(dir_name))
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def read_mol(input_mol, folder_path, is_core=False):
    """
    Checks by which means the input is provided and converts it into a molecule.
    Returns a list of PLAMS molecules.
    """
    # Creates a dictionary of file extensions
    extension_dict = {'xyz': read_mol_xyz, 'pdb': read_mol_pdb, 'mol': read_mol_mol,
                      'smiles': read_mol_smiles, 'folder': read_mol_folder, 'txt': read_mol_txt,
                      'xlsx': read_mol_excel, 'plams_mol': read_mol_plams, 'rdmol': read_mol_rdkit}

    # Reads the input molecule(s), the method depending on the nature of the file extension
    # Returns a list of dictionaries
    input_mol = [read_mol_extension(mol, folder_path, is_core) for mol in input_mol]
    input_mol = dict_concatenate(input_mol)

    # Turns input_mol into a list if it is a string
    if isinstance(input_mol, str):
        input_mol = [input_mol]

    # Create a list of PLAMS molecules
    mol_list = []
    for mol in input_mol:
        try:
            read = extension_dict[input_mol[mol][0]]
            mol_list.append(read(mol, input_mol[mol][1]))
        except KeyError as ex:
            print(str(type(ex).__name__) + ':\t' + str(ex) + '\n')
    mol_list = list(itertools.chain(*mol_list))

    # Raises an error if mol_list is empty
    if mol_list:
        return mol_list
    else:
        core_ligand = {True: 'cores', False: 'ligands'}
        raise IndexError('No valid input ' + core_ligand[is_core] + ' were found, aborting run')


def read_mol_extension(mol_name, folder_path, is_core=False):
    """
    Identifies the filetypes used in the input molecules.
    Returns the input molecule (key), the filetype (arg1) and optional arguments (arg2)
    """
    # Check if mol_name contains optional keyword arguments
    if isinstance(mol_name, list) or isinstance(mol_name, tuple):
        kwarg = mol_name[1]
        mol_name = mol_name[0]
    else:
        kwarg = {}
    kwarg.update({'folder_path': folder_path, 'is_core': is_core})

    # Identify the filetype of mol_name
    if mol_name is None:
        mol_name = ''
    mol_path = os.path.join(folder_path, mol_name)
    if os.path.isfile(mol_path):
        return {mol_name: [mol_name.rsplit('.', 1)[-1], kwarg]}
    elif os.path.isdir(mol_path):
        return {mol_name: ['folder', kwarg]}
    elif isinstance(mol_name, Molecule):
        return {mol_name: ['plams_mol', kwarg]}
    elif isinstance(mol_name, Chem.rdchem.Mol):
        return {mol_name: ['rdmol', kwarg]}
    else:
        return {mol_name: ['smiles', kwarg]}


def read_mol_xyz(mol_name, kwarg):
    """
    Read an .xyz file
    """
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        mol = Molecule(mol_path, inputformat='xyz')
        mol_name = mol_name.rsplit('.', 1)[0]
        if kwarg.get('guess_bonds') is None or kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_xyz.__code__, ex, mol_path, [mol_name, kwarg])


def read_mol_pdb(mol_name, kwarg):
    """
    Read a .pdb file
    """
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        mol = molkit.readpdb(mol_path)
        mol_name = mol_name.rsplit('.', 1)[0]
        if kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_pdb.__code__, ex, mol_path, [mol_name, kwarg])


def read_mol_mol(mol_name, kwarg):
    """
    Read a .mol file
    """
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        mol = molkit.from_rdmol(Chem.MolFromMolFile(mol_path, removeHs=False))
        mol_name = mol_name.rsplit('.', 1)[0]
        if kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_mol.__code__, ex, mol_path, [mol_name, kwarg])


def read_mol_smiles(mol_name, kwarg):
    """
    Read a SMILES string
    """
    try:
        mol = molkit.from_smiles(mol_name)
        mol_name = mol_name.replace('/', '').replace('\\', '')
        if kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_smiles.__code__, ex, mol_name, [mol_name, kwarg])


def read_mol_plams(mol_name, kwarg):
    """
    Read a PLAMS molecule
    """
    try:
        mol = mol_name
        if mol.properties.name:
            mol_name = mol.properties.name
        else:
            mol_name = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))
        if kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_plams.__code__, ex, mol_name, [mol_name, kwarg])


def read_mol_rdkit(mol_name, kwarg):
    """
    Read a RDKit molecule
    """
    try:
        mol = molkit.from_rdmol(mol_name)
        if mol.properties.name:
            mol_name = mol.properties.name
        else:
            mol_name = Chem.MolToSmiles(Chem.RemoveHs(mol_name))
        if kwarg.get('guess_bonds'):
            mol.guess_bonds()
        set_prop(mol, mol_name, kwarg['folder_path'], kwarg['is_core'])
        return [mol]
    except (Exception) as ex:
        return print_exception(read_mol_rdkit.__code__, ex, mol_name, [mol_name, kwarg])


def read_mol_folder(mol_name, kwarg):
    """
    Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder
    """
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        file_list = [[file, kwarg] for file in os.listdir(mol_path)]
        return read_mol(file_list, kwarg['folder_path'], kwarg['is_core'])
    except (Exception, errors.PlamsError) as ex:
        return print_exception(read_mol_folder.__code__, ex, mol_path, [mol_name, kwarg])


def read_mol_txt(mol_name, kwarg):
    """
    Read a plain text file containing one or more SMILES strings
    """
    kwarg = dict_defaults(kwarg, {'row': 0, 'column': 0})
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        with open(mol_path, 'r') as file:
            mol_list = file.read().splitlines()
        mol_list = [mol.split()[kwarg['column']] for mol in mol_list[kwarg['row']:] if mol]
        mol_list = [[mol, kwarg] for mol in mol_list if len(mol) >= 2]
        return read_mol(mol_list, kwarg['folder_path'], kwarg['is_core'])
    except (Exception) as ex:
        return print_exception(read_mol_txt.__code__, ex, mol_path, [mol_name, kwarg])


def read_mol_excel(mol_name, kwarg):
    """
    Read a plain text file containing one or more SMILES strings
    """
    kwarg = dict_defaults(kwarg, {'row': 0, 'column': 0, 'sheet_name': 'Sheet1'})
    mol_path = os.path.join(kwarg['folder_path'], mol_name)
    try:
        mol_list = pd.read_excel(mol_path, sheet_name=kwarg['sheet_name'])
        mol_list = [[mol, kwarg] for mol in mol_list[kwarg['column']][kwarg['row']:]]
        return read_mol(mol_list, kwarg['folder_path'], kwarg['is_core'])
    except (Exception) as ex:
        return print_exception(read_mol_excel.__code__, ex, mol_path, [mol_name, kwarg])


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
