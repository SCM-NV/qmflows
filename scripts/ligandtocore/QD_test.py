import QD_import_export as QD_inout


def test_input():
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

    assert [formula_list[1] == formula for formula in formula_list[2:]]
    assert [distance_list[1] == distance for distance in distance_list[2:]]
    assert [angle_list[1] == angle for angle in angle_list[2:]]
