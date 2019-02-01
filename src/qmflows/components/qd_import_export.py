__all__ = ['read_mol', 'set_prop', 'create_dir']

import os
import itertools
import pandas as pd

from scm.plams import Molecule
from scm.plams.core import errors
import scm.plams.interfaces.molecule.rdkit as molkit

from rdkit import Chem

from .qd_functions import get_time


def read_mol(input_mol, path, is_core=False):
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
    input_mol = [read_mol_defaults(mol, path, is_core) for mol in input_mol]
    input_mol = [read_mol_extension(mol) for mol in input_mol]
    input_mol = dict_concatenate(input_mol)

    # Create a list of PLAMS molecules
    mol_list = []
    for mol in input_mol:
        mol_dict = input_mol[mol]

        # Check if the provided key is available in extension_dict
        try:
            read = extension_dict[mol_dict['file_type']]
        except KeyError as ex:
            print(get_time() + str(type(ex).__name__) + ':\t' + str(ex) + '\n')

        # Convert mol into either a list or PLAMS molecule
        if read:
            mol = read(mol, mol_dict)
            # if mol is a PLAMS molecule
            if isinstance(mol, Molecule):
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
            elif isinstance(mol, list):
                mol_list += mol

    return mol_list


def read_mol_defaults(mol, path, is_core=False):
    """
    Assign all default and/or user-specified arguments to mol.

    mol <str>, <none> or <list>[<str>, <dict>]: The input molecule.
    path <str>: The path to mol.
    is_core <bool>: If the input mol is a core (True) or ligand (False).

    return <dict>: A dictionary containing mol (key) and all default and/or user-specified args.
    """
    # A dictionary of default arguments
    kwarg = {'file_type': '',
             'name': '',
             'row': 0,
             'column': 0,
             'sheet_name': 'Sheet1',
             'guess_bonds': False,
             'is_core': is_core,
             'path': path,
             'user_args': mol,
             'core_indices': '',
             'ligand_indices': ''}

    # Update the dictionary of default arguments based on used-specified arguments
    if isinstance(mol, (list, tuple)):
        for item in mol[1:]:
            kwarg.update(item)
        mol = mol[0]
    if mol is None or mol == 'None':
        mol = ''
    kwarg['name'] = mol

    return {mol: kwarg}


def read_mol_extension(mol_dict):
    """
    Identifies the filetype of mol and updates the dictionary accordingly.

    mol_dict <dict>: A dictionary containing mol (key) and all default and/or user-specified args.

    return <dict>: A dictionary containing mol (key) and all default and/or user-specified args,
        now with an updated name and file_type.
    """
    mol = list(mol_dict.keys())[0]
    try:
        mol_path = os.path.join(mol_dict[mol]['path'], mol_dict[mol]['name'])
        if os.path.isfile(mol_path):
            file_type = mol.rsplit('.', 1)[-1]
            name = mol.rsplit('.', 1)[0]
        elif os.path.isdir(mol_path):
            file_type = 'folder'
            name = mol
        elif os.path.isfile(mol):
            file_type = mol.rsplit('.', 1)[-1]
            name = mol.rsplit('.', 1)[0].rsplit('/', 1)[-1].rsplit('\\', 1)[-1]
        elif os.path.isdir(mol):
            file_type = 'folder'
            name = mol
        else:
            file_type = 'smiles'
            name = smiles_name(mol)
    except TypeError:
        if isinstance(mol, Molecule):
            file_type = 'plams_mol'
            if not mol.properties.name:
                name = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))
            else:
                name = mol.properties.name
        elif isinstance(mol, Chem.rdchem.Mol):
            file_type = 'rdmol'
            name = Chem.MolToSmiles(Chem.RemoveHs(mol))

    # Update mol_dict
    mol_dict[mol]['file_type'] = file_type
    mol_dict[mol]['name'] = name

    return mol_dict


def smiles_name(smiles_str):
    """
    Turn a SMILES string into an acceptable filename.

    smiles_str <str>: A SMILES string.

    return <str>: A filename based on a SMILES string.
    """
    name = smiles_str.replace('(', '[').replace(')', ']')
    cis_trans = [item for item in smiles_str if item == '/' or item == '\\']
    if cis_trans:
        cis_trans = [item + cis_trans[i*2+1] for i, item in enumerate(cis_trans[::2])]
        cis_trans_dict = {'//': 'trans-', '/\\': 'cis-'}
        for item in cis_trans[::-1]:
            name = cis_trans_dict[item] + name
        name = name.replace('/', '').replace('\\', '')

    return name


def read_mol_xyz(mol, mol_dict):
    """
    Read an .xyz file
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'] + '.xyz')
        return Molecule(mol_path, inputformat='xyz')
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_xyz.__code__, ex, mol_path)


def read_mol_pdb(mol, mol_dict):
    """
    Read a .pdb file
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'] + '.pdb')
        return molkit.readpdb(mol_path)
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_pdb.__code__, ex, mol_path)


def read_mol_mol(mol, mol_dict):
    """
    Read a .mol file
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'] + '.mol')
        return molkit.from_rdmol(Chem.MolFromMolFile(mol_path, removeHs=False))
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_mol.__code__, ex, mol_path)


def read_mol_smiles(mol, mol_dict):
    """
    Read a SMILES string
    """
    try:
        return molkit.from_smiles(mol)
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_smiles.__code__, ex, mol)


def read_mol_plams(mol, mol_dict):
    """
    Read a PLAMS molecule
    """
    try:
        return mol
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_plams.__code__, ex, mol)


def read_mol_rdkit(mol, mol_dict):
    """
    Read a RDKit molecule
    """
    try:
        return molkit.from_rdmol(mol)
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_rdkit.__code__, ex, mol)


def read_mol_folder(mol, mol_dict):
    """
    Read all files (.xyz, .pdb, .mol, .txt or further subfolders) within a folder
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'])
        file_list = [[file, mol_dict] for file in os.listdir(mol_path)]
        return read_mol(file_list, mol_dict['path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_folder.__code__, ex, mol_path)


def read_mol_txt(mol, mol_dict):
    """
    Read a plain text file containing one or more SMILES strings
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'] + '.txt')
        with open(mol_path, 'r') as file:
            file_list = file.read().splitlines()
        file_list = [file.split()[mol_dict['column']] for file in file_list[mol_dict['row']:] if
                     file]
        file_list = [[file, mol_dict] for file in file_list if len(file) >= 2]
        return read_mol(file_list, mol_dict['path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_txt.__code__, ex, mol_path)


def read_mol_excel(mol, mol_dict):
    """
    Read a plain text file containing one or more SMILES strings
    """
    try:
        mol_path = os.path.join(mol_dict['path'], mol_dict['name'] + '.xlsx')
        file_list = list(pd.read_excel(mol_path, sheet_name=mol_dict['sheet_name']))
        file_list = [[file, mol_dict] for file in file_list[mol_dict['column']][mol_dict['row']:]]
        return read_mol(file_list, mol_dict['path'], mol_dict['is_core'])
    except (Exception, errors.PlamsError) as ex:
        print_exception(read_mol_excel.__code__, ex, mol_path)


def set_prop(mol, mol_dict):
    """
    Set molecular and atomic properties
    """
    if mol_dict.get('is_core'):
        residue_name = 'COR'
        mol.properties.name = 'Core_' + mol_dict['name']
        mol.properties.dummies = mol_dict['core_indices']
    else:
        residue_name = 'LIG'
        mol.properties.name = 'Ligand_' + mol_dict['name']
        mol.properties.dummies = mol_dict['ligand_indices']

    mol.properties.path = mol_dict['path']

    # Prepare a list of letters for pdb_info.Name
    alphabet = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet = [i + j for i in alphabet for j in alphabet]

    # Create a dictionary of elements and their formal atomic charge
    elements_dict = set_prop_dict()

    # Set the atomic properties
    for i, atom in enumerate(mol):
        set_prop_atom(atom, alphabet[i], residue_name, elements_dict)

    if not mol.properties.smiles:
        mol.properties.smiles = Chem.MolToSmiles(Chem.RemoveHs(molkit.to_rdmol(mol)))


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
            total_bonds = int(sum([bond.order for bond in atom.bonds]))
            default_charge = elements_dict[atom.symbol]
            sign = int(-1 * default_charge / abs(default_charge))
            atom.properties.charge = default_charge + sign*total_bonds

            # Update formal atomic charges for hypervalent atoms
            if total_bonds > abs(default_charge):
                if total_bonds is abs(default_charge) + 2:
                    atom.properties.charge += sign*2
                elif total_bonds is abs(default_charge) + 4:
                    atom.properties.charge += sign*4
                elif total_bonds >= abs(default_charge) + 6:
                    atom.properties.charge += sign*6
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
    name = mol.properties.name
    mol_path = os.path.join(mol.properties.path, name)
    molkit.writepdb(mol, mol_path + '.pdb')
    mol.write(mol_path + '.xyz')
    print(get_time() + str(message) + str(name) + '.pdb')


def print_exception(func, ex, name):
    """
    Manages the printing of exceptions upon failing to import a molecule.
    """
    extension_dict = {'read_mol_xyz': '.xyz file', 'read_mol_pdb': '.pdb file',
                      'read_mol_mol': '.mol file', 'read_mol_smiles': 'SMILES string',
                      'read_mol_folder': 'folder', 'read_mol_txt': '.txt file',
                      'read_mol_excel': '.xlsx file', 'read_mol_plams': 'PLAMS molecule',
                      'read_mol_rdkit': 'RDKit molecule'}
    print(get_time() + str(type(ex).__name__), str(ex))
    print(get_time() + 'function:', str(func.co_name) + str(func.co_varnames[:func.co_argcount]))
    # print('\t\targ1(mol):', type(mol), mol)
    # print('\t\targ2(mol_dict):', type(mol_dict), mol_dict)
    print(get_time() + 'Warning:', name, 'not recognized as a valid',
          extension_dict[func.co_name], '\n')
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
