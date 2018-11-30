import os
import pytest
import yaml

from scm.plams import Molecule
import qmflows.qd as QD
from qmflows.components.qd_import_export import (read_mol, read_mol_xyz, read_mol_pdb, read_mol_mol)

path = os.path.abspath('test/test_QD/test_QD_files')
# path = os.path.join(os.getcwd(), 'test_qd_files')


def test_read_mol_1():
    """
    Test if the AcOH.xyz, .pdb and .mol files return molecules with molecular formulas,
    distances and angles.
    i.e. are all internal coordinates identical?
    """
    xyz = read_mol_xyz('AcOH.xyz', {'path': path, 'name': 'AcOH'})
    pdb = read_mol_pdb('AcOH.pdb', {'path': path, 'name': 'AcOH'})
    mol = read_mol_mol('AcOH.mol', {'path': path, 'name': 'AcOH'})
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

    assert [formula_list[0] is formula for formula in formula_list[1:]]
    assert [distance_list[0] is distance for distance in distance_list[1:]]
    assert [angle_list[0] is angle for angle in angle_list[1:]]


def test_read_mol_2():
    """
    Check if 5 valid SMILES strings return 5 plams molecules
    """
    input_mol = ['OC', 'OCC', 'OCCC', 'OCCCC', 'OCCCCC']
    mol_list = read_mol(input_mol, path)

    assert isinstance(mol_list, list)
    assert len(mol_list) is len(input_mol)
    for mol in mol_list:
        assert isinstance(mol, Molecule)


def test_input():
    path = os.path.abspath('test/test_QD')
    # path = os.getcwd()
    input_cores = yaml.load("""
    -   - Cd68Se55.xyz
        - guess_bonds: False
    """)

    input_ligands = yaml.load("""
    - AcOH.mol
    - AcOH.pdb
    - AcOH.txt
    - AcOH.xyz
    - O=C(O)C
    """)

    argument_dict = yaml.load("""
    dir_name_list: [test_qd_files, test_qd_files, test_qd_files]
    dummy: Cl
    database_name: [ligand_database.xlsx, QD_database.xlsx]
    use_database: True
    core_opt: False
    ligand_opt: True
    ligand_crs: False
    qd_opt: False
    qd_int: True
    maxiter: 500
    split: True
    """)

    qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
    QD.prep(input_ligands, input_cores, path, argument_dict)
    formula_set = set([qd.get_formula() for qd in qd_list])

    assert isinstance(qd_list, list)
    assert len(qd_list) is 5
    assert len(formula_set) is 1
    for qd in qd_list:
        assert isinstance(qd, Molecule)

    print(True)
    argument_dict['use_database'] = False
    QD.prep(input_ligands, input_cores, path, argument_dict)
    argument_dict['ligand_opt'] = False
    QD.prep(input_ligands, input_cores, path, argument_dict)
    argument_dict['split'] = False
    QD.prep(input_ligands, input_cores, path, argument_dict)
    argument_dict['dummy'] = 'Cd'
    QD.prep(input_ligands, input_cores, path, argument_dict)
    argument_dict['qd_int'] = False
    QD.prep(input_ligands, input_cores, path, argument_dict)

    os.remove(os.path.join(path, 'Ligand_database.xlsx'))
    exclusion = ['AcOH.mol', 'AcOH.pdb', 'AcOH.txt', 'AcOH.xyz', 'Cd68Se55.xyz']
    path = os.path.join(path, 'test_qd_files')
    for file in reversed(os.listdir(path)):
        if file not in exclusion:
            os.remove(os.path.join(path, file))
