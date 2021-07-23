"""Test reactivity utilities."""

import pytest
from qmflows import Settings
from scm.plams import Molecule
from qmflows.test_utils import PATH_MOLECULES, HAS_RDKIT
from assertionlib import assertion

if HAS_RDKIT:
    from qmflows.components import reactivity

mol = Molecule(PATH_MOLECULES / "ethylene.xyz")


@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_Distance():
    """Test Distance reactivity functionality."""
    d = reactivity.Distance(1, 2)
    s = d.get_settings(mol=mol)

    assertion.isinstance(s, Settings)


@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_Angle():
    """Test Angle reactivity functionality."""
    ang = reactivity.Angle(1, 2, 3)
    s = ang.get_settings(mol=mol)

    assertion.isinstance(s, Settings)


@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_Dihedral():
    """Test Angle reactivity functionality."""
    dihed = reactivity.Dihedral(1, 2, 3, 4)
    s = dihed.get_settings(mol=mol)

    assertion.isinstance(s, Settings)
