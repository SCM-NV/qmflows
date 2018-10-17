import itertools
import os
import pytest

from scm.plams import Molecule
import QD_import_export as QD_inout


def test_read_mol_1():
    """
    Test if the AcOH.xyz, .pdb and .mol files return molecules with molecular formulas,
    distances and angles.
    i.e. are all internal coordinates identical?
    """
    mol_name = 'AcOH'
    kwarg = {'is_core': False, 'folder_path': '/Users/basvanbeek/Documents/GitHub/' +
                                              'qmflows/scripts/ligandtocore/test_mol'}

    xyz = QD_inout.read_mol_xyz(mol_name + '.xyz', kwarg)[0]
    pdb = QD_inout.read_mol_pdb(mol_name + '.pdb', kwarg)[0]
    mol = QD_inout.read_mol_mol(mol_name + '.mol', kwarg)[0]
    mol_list = [xyz, pdb, mol]

    atom_list = [[[at1, at2] for at1 in mol for at2 in mol if at1 != mol[1] and at2 != mol[1]] for
                 mol in mol_list]
    formula_list = [mol.get_formula() for mol in mol_list]
    distance_list = [[mol[1].distance_to(atom) for atom in mol if not mol[1]] for mol in mol_list]
    angle_list = []
    for i, mol in enumerate(mol_list):
        for atom in atom_list[i]:
            if atom[0] != atom[1]:
                angle_list.append([mol[1].angle(atom[0], atom[1], result_unit='degree')])

    assert isinstance(mol_list, list)
    assert len(mol_list) == 3
    for mol in mol_list:
        assert isinstance(mol, Molecule)
        assert len(mol.properties) >= 4

    assert [formula_list[0] is formula for formula in formula_list[1:]]
    assert [distance_list[0] is distance for distance in distance_list[1:]]
    assert [angle_list[0] is angle for angle in angle_list[1:]]


def test_read_mol_2():
    """
    Checks parsing if a file with a given extensions returns:
        1. a list
        2. if the list has 1 entry
        3. if that entry is a plams molecule
    """
    mol_name = 'AcOH'
    kwarg = {'is_core': False, 'folder_path': '/Users/basvanbeek/Documents/GitHub/' +
                                              'qmflows/scripts/ligandtocore/test_mol'}
    extension_dict = {'xyz': QD_inout.read_mol_xyz, 'pdb': QD_inout.read_mol_pdb,
                      'mol': QD_inout.read_mol_mol, 'smiles': QD_inout.read_mol_smiles,
                      'folder': QD_inout.read_mol_folder, 'txt': QD_inout.read_mol_txt,
                      'xlsx': QD_inout.read_mol_excel, 'plams_mol': QD_inout.read_mol_plams,
                      'rdmol': QD_inout.read_mol_rdkit}
    key_list = list(extension_dict.keys())
    for item in key_list:
        path = os.path.join(kwarg['folder_path'], mol_name + '.' + item)
        if os.path.exists(path):
            xyz = [extension_dict[key](mol_name + '.' + item, kwarg) for key in key_list]
            xyz = list(itertools.chain(*xyz))

            assert isinstance(xyz, list)
            assert len(xyz) == 1
            assert isinstance(xyz[0], Molecule)
            assert len(xyz[0].properties) >= 4


def test_read_mol_3():
    """
    Check if 5 valid SMILES strings return 5 plams molecules
    """
    input_mol = ['OC', 'OCC', 'OCCC', 'OCCCC', 'OCCCCC']
    folder_path = '/Users/basvanbeek/Documents/GitHub/qmflows/scripts/ligandtocore/test_mol'
    mol_list = QD_inout.read_mol(input_mol, folder_path=folder_path)

    assert isinstance(mol_list, list)
    assert len(mol_list) is len(input_mol)
    for mol in mol_list:
        assert isinstance(mol, Molecule)
        assert len(mol.properties) >= 4


def test_read_mol_4():
    """
    Check if 2 invalid SMILES strings return an IndexError
    """
    input_mol = ['dwefwefqe', 'fqwdwq']
    folder_path = '/Users/basvanbeek/Documents/GitHub/qmflows/scripts/ligandtocore/test_mol'
    with pytest.raises(IndexError):
        assert QD_inout.read_mol(input_mol, folder_path=folder_path)
