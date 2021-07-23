"""Test molecule serialization."""
from qmflows import packages
from qmflows.test_utils import PATH
from scm.plams import Molecule

WATER = Molecule(PATH / "water.xyz")
WATER.guess_bonds()


def test_SerMolecule():
    """Test molecule serialization."""
    mol = WATER
    registry = packages.registry()
    encoded_molecule = registry.deep_encode(mol)
    decoded_molecule = registry.deep_decode(encoded_molecule)
    assert len(mol) == len(decoded_molecule)
