"""Test input keyword with special meaning."""

from os.path import basename

import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, adf, dftb, orca
from qmflows.fileFunctions import yaml2Settings
from qmflows.packages.cp2k_package import CP2K
from qmflows.packages.cp2k_mm import CP2KMM
from qmflows.test_utils import PATH_MOLECULES, PATH, get_mm_settings

indices = [3, 4, 5, 6]
adf_const = {
    'adf': {'constraints': {'atom 1': '', 'atom 2': ''},
            'geometry': {'optim': 'cartesian'}}}
dftb_const = {
    'dftb': {'constraints': {'atom 1': '', 'atom 2': ''},
             'geometry': {'optim': 'cartesian'}}}
orca_const = {
    'orca': {'geom': {'Constraints': {'_end': '{ C 0 C }{ C 1 C }'}}}}


def test_freeze_with_adf():
    """Test the freeze keyword with ADF."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ["C", "O"], 'specific': adf_const})
    assertion.eq(str(adf.generic2specific(s, mol)), str(expected_settings))
    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': adf_const})
    assertion.eq(str(adf.generic2specific(s, mol)), str(expected_settings))


def test_selected_atoms_with_adf():
    """Test the ADF atoms selection features."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': adf_const})

    assertion.eq(str(adf.generic2specific(s, mol)), str(expected_settings))
    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': adf_const})
    assertion.eq(str(adf.generic2specific(s, mol)), str(expected_settings))


def test_freeze_with_dftb():
    """Test freeze keyword for DFTB."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ['C', 'O'], 'specific': dftb_const})
    assertion.eq(str(dftb.generic2specific(s, mol)), str(expected_settings))
    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': dftb_const})
    assertion.eq(str(dftb.generic2specific(s, mol)), str(expected_settings))


def test_selected_atoms_with_dftb():
    """Test the DFTB selection atoms funcionality."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': dftb_const})
    assertion.eq(dftb.generic2specific(s, mol), expected_settings)
    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': dftb_const})
    assertion.eq(str(dftb.generic2specific(s, mol)), str(expected_settings))


def test_freeze_with_orca():
    """Test freeze keyword for orca."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.freeze = ["C", "O"]
    expected_settings = Settings(
        {'freeze': ['C', 'O'], 'specific': orca_const})
    assertion.eq(str(orca.generic2specific(s, mol)), str(expected_settings))

    s.freeze = [1, 2]
    expected_settings = Settings(
        {'freeze': [1, 2], 'specific': orca_const})

    assertion.eq(str(orca.generic2specific(s, mol)), str(expected_settings))


def test_selected_atoms_with_orca():
    """Test atom selection features for Orca."""
    mol = molkit.from_smiles('CO')
    s = Settings()
    s.selected_atoms = ["H"]
    expected_settings = Settings(
        {'selected_atoms': ['H'], 'specific': orca_const})
    assertion.eq(str(orca.generic2specific(s, mol)), str(expected_settings))

    s.selected_atoms = indices
    expected_settings = Settings(
        {'selected_atoms': indices, 'specific': orca_const})

    assertion.eq(str(orca.generic2specific(s, mol)), str(expected_settings))


def test_cp2k_special_keywords():
    """Test the translation from settings to CP2K specific keywords."""
    abc = [5.958, 7.596, 15.610]
    angles = [81.25, 86.56, 89.80]
    ETHYLENE = Molecule(PATH_MOLECULES / "ethylene.xyz")

    s = Settings()
    s.cell_parameters = abc
    s.cell_angles = angles
    s.periodic = None

    # apply transformations
    CP2K.handle_special_keywords(s, "cell_parameters", abc, ETHYLENE)
    CP2K.handle_special_keywords(s, "cell_angles", angles, ETHYLENE)
    CP2K.handle_special_keywords(s, "periodic", None, ETHYLENE)

    # compare with the reference
    ref = Settings()
    ref.specific.cp2k.force_eval.subsys.cell.ABC = " [angstrom] 5.958 7.596 15.61"
    ref.specific.cp2k.force_eval.subsys.cell.ALPHA_BETA_GAMMA = "81.25 86.56 89.8"
    ref.specific.cp2k.force_eval.subsys.cell.periodic = None
    assertion.eq(s.specific, ref.specific)


def test_cp2k_mm_keywords():
    """Test the translation from settings to CP2KMM specific keywords."""
    s = get_mm_settings()

    CP2KMM.prerun(None, s, None)
    CP2KMM.handle_special_keywords(s, 'psf', s.psf, None)
    CP2KMM.handle_special_keywords(s, 'prm', s.prm, None)
    CP2KMM.handle_special_keywords(s, 'lennard_jones', s.lennard_jones, None)
    CP2KMM.handle_special_keywords(s, 'charge', s.charge, None)
    CP2KMM.handle_special_keywords(s, 'periodic', s.periodic, None)

    # Change absolute to relative path names for the purpose of testing
    s.specific.cp2k.force_eval.subsys.topology.conn_file_name = basename(s.specific.cp2k.force_eval.subsys.topology.conn_file_name)
    s.specific.cp2k.force_eval.mm.forcefield.parm_file_name = basename(s.specific.cp2k.force_eval.mm.forcefield.parm_file_name)

    with open(PATH / 'cp2k_mm_special_keyword.yaml', 'rb') as f:
        ref = yaml2Settings(f)
    assertion.eq(s.specific, ref)
