from qmflows import Settings, gamess, adf, dftb, orca, molkit

def test_freeze_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings({'freeze': ['C', 'O'], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.freeze = [1,2]
    expected_settings = Settings({'freeze': [1, 2], 'specific': {'gamess': {'statpt': 'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)

def test_selected_atoms_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings({'selected_atoms': ['H'], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.selected_atoms = [3,4,5,6]
    expected_settings = Settings({'selected_atoms': [3,4,5,6], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)

def test_freeze_with_adf():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings({'freeze': ["C", "O"], 'specific':
                                    {'adf': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)
    s.freeze = [1, 2]
    expected_settings = Settings({'freeze': [1, 2], 'specific':
                                    {'adf': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)

def test_selected_atoms_with_adf():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings({'selected_atoms': ['H'], 'specific':
                                    {'adf': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)
    s.selected_atoms = [3,4,5,6]
    expected_settings = Settings({'selected_atoms': [3,4,5,6], 'specific':
                                    {'adf': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)

def test_freeze_with_dftb():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings({'freeze': ['C', 'O'], 'specific':
                                    {'dftb': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)
    s.freeze = [1, 2]
    expected_settings = Settings({'freeze': [1, 2], 'specific':
                                    {'dftb': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)

def test_selected_atoms_with_dftb():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings({'selected_atoms': ['H'], 'specific':
                                    {'dftb': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert dftb.generic2specific(s, mol) == expected_settings
    s.selected_atoms = [3,4,5,6]
    expected_settings = Settings({'selected_atoms': [3,4,5,6], 'specific':
                                    {'dftb': {'constraints': {'atom 1': '', 'atom 2': ''}, 'geometry': {'optim': 'cartesian'}}}})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)

def test_freeze_with_orca():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings({'freeze': ['C', 'O'], 'specific':
        {'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)
    s.freeze = [1, 2]
    expected_settings = Settings({'freeze': [1, 2], 'specific':
        {'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)

def test_selected_atoms_with_orca():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings({'selected_atoms': ['H'], 'specific':
        {'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)
    s.selected_atoms = [3,4,5,6]
    expected_settings = Settings({'selected_atoms': [3,4,5,6], 'specific':
        {'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)
