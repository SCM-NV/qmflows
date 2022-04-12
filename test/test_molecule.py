"""Test molecule serialization."""

from qmflows.packages._packages import registry
from qmflows.test_utils import PATH
from scm.plams import Molecule

WATER = Molecule(PATH / "water.xyz")
WATER.guess_bonds()


def test_SerMolecule():
    """Test molecule serialization."""
    mol = WATER
    reg = registry()
    encoded_molecule = reg.deep_encode(mol)
    decoded_molecule = reg.deep_decode(encoded_molecule)
    assert len(mol) == len(decoded_molecule)
