"""Test input keyword with special meaning."""
import numpy as np
import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, adf, dftb, orca
from qmflows.packages.cp2k_package import CP2K
from qmflows.packages.orca import ORCA
from qmflows.packages.SCM import ADF
from qmflows.parsers.adf_parser import kfreader
from qmflows.parsers.orca_parser import parse_hessian
from qmflows.test_utils import PATH, PATH_MOLECULES

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


def test_orca_init_hessian():
    """Test the translation from settings to CP2K specific keywords."""
    # Read Hessian from DFTB
    PATH_RKF = PATH / "output_dftb" / "dftb_freq" / "dftb.rkf"
    assertion.truth(PATH_RKF.exists())
    hess = kfreader(PATH_RKF, section="AMSResults", prop="Hessian")
    water = molkit.from_smiles('[OH2]', forcefield='mmff')

    # Tess Hessian initialization
    s = Settings()
    hess = np.array(hess).reshape(9, 9)
    s.inithess = hess
    ORCA.handle_special_keywords(s, "inithess", hess, water)

    # Test that the hessian is readable
    new_hess = parse_hessian(s.specific.orca.geom.InHessName)
    new_hess = np.array(new_hess, dtype=np.float)
    assertion.truth(np.allclose(hess, new_hess, atol=1e-5))


def test_orca_constrains():
    """Test geometry constrains in orca."""
    ethylene = PATH_MOLECULES / "molecules" / "ethylene.xyz"

    # Test distance constrains
    s = Settings()
    s.constraint['dist 1 2'] = 1.1  # Constrain C-H bond to 1.1 Angstrom
    ORCA.handle_special_keywords(
        s, "constraint", Settings({'dist 1 2': 1.1}), ethylene)
    assertion.eq(s.specific.orca.geom.Constraints._end, "{ B 0 1 1.10 C }")

    # Test angle constrains
    s = Settings()
    s.constraint['angle 1 2 3'] = 109.5  # Constrain H-C-H to 109.5
    ORCA.handle_special_keywords(
        s, "constraint", Settings({'angle 1 2 3': 109.5}), ethylene)
    assertion.eq(s.specific.orca.geom.Constraints._end, "{ A 0 1 2 109.50 C }")

    s = Settings()
    s.constraint['dihed 1 2 3 4'] = 180  # Constrain Dihedral to 180
    ORCA.handle_special_keywords(
        s, "constraint", Settings({'dihed 1 2 3 4': 180}), ethylene)
    assertion.eq(s.specific.orca.geom.Constraints._end,
                 "{ D 0 1 2 3 180.00 C }")


def test_adf_constrains():
    """Test the geometry constrains in ADF."""
    ethylene = PATH_MOLECULES / "molecules" / "ethylene.xyz"

    # Test distance constrains
    s = Settings()
    s.constraint['dist 1 2'] = 1.1  # Constrain C-H bond to 1.1 Angstrom
    ADF.handle_special_keywords(
        s, "constraint", Settings({'dist 1 2': 1.1}), ethylene)
    assertion.eq(s.specific.adf.constraints, Settings({"dist 1 2": 1.1}))

    # Test distance angles
    s = Settings()
    s.constraint['angle 1 2 3'] = 109.5  # Constrain C-H bond to 1.1 Angstrom
    ADF.handle_special_keywords(
        s, "constraint", Settings({'angle 1 2 3': 109.5}), ethylene)
    assertion.eq(s.specific.adf.constraints, Settings({'angle 1 2 3': 109.5}))

    # Test distance dihedral
    s = Settings()
    s.constraint['dihe 1 2 3 4'] = 180  # Constrain C-H bond to 1.1 Angstrom
    ADF.handle_special_keywords(
        s, "constraint", Settings({'dihed 1 2 3 4': 180}), ethylene)
    assertion.eq(s.specific.adf.constraints, Settings({'dihed 1 2 3 4': 180}))
