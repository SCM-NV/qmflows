from qmworks import Settings, gamess, molkit

def test_freeze_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings({'freeze': ['C', 'O'], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.freeze = [0,1]
    expected_settings = Settings({'freeze': [0, 1], 'specific': {'gamess': {'statpt': 'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)

def test_selected_atoms_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings({'selected_atoms': ['H'], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.selected_atoms = [2,3,4,5]
    expected_settings = Settings({'selected_atoms': [2,3,4,5], 'specific': {'gamess': {'statpt':'IFREEZ(1)=1,2,3,4,5,6'}}})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)

