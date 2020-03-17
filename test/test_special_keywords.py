"""Test input keyword with special meaning."""
from qmflows import (Settings, gamess, adf, dftb, orca)
import scm.plams.interfaces.molecule.rdkit as molkit

indices = [3, 4, 5, 6]
gamess_sett = {'gamess': {'statpt': 'IFREEZ(1)=1,2,3,4,5,6'}}
adf_const = {
    'adf': {'constraints': {'atom 1': '', 'atom 2': ''},
            'geometry': {'optim': 'cartesian'}}}
dftb_const = {
    'dftb': {'constraints': {'atom 1': '', 'atom 2': ''},
             'geometry': {'optim': 'cartesian'}}}
orca_const = {
    'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}


def test_freeze_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ['C', 'O'], 'specific': gamess_sett})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': gamess_sett})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)


def test_selected_atoms_with_gamess():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': gamess_sett})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)
    s = Settings()
    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': gamess_sett})
    assert str(gamess.generic2specific(s, mol)) == str(expected_settings)


def test_freeze_with_adf():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ["C", "O"], 'specific': adf_const})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)
    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': adf_const})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)


def test_selected_atoms_with_adf():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': adf_const})

    assert str(adf.generic2specific(s, mol)) == str(expected_settings)
    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': adf_const})
    assert str(adf.generic2specific(s, mol)) == str(expected_settings)


def test_freeze_with_dftb():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ['C', 'O'], 'specific': dftb_const})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)
    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': dftb_const})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)


def test_selected_atoms_with_dftb():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': dftb_const})
    assert dftb.generic2specific(s, mol) == expected_settings
    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': dftb_const})
    assert str(dftb.generic2specific(s, mol)) == str(expected_settings)


def test_freeze_with_orca():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ['C', 'O'], 'specific': orca_const})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)

    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': orca_const})

    assert str(orca.generic2specific(s, mol)) == str(expected_settings)


def test_selected_atoms_with_orca():
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': orca_const})
    assert str(orca.generic2specific(s, mol)) == str(expected_settings)

    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': orca_const})

    assert str(orca.generic2specific(s, mol)) == str(expected_settings)
