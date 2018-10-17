from qmflows import Settings
from qmflows.components import reactivity
from scm.plams import Molecule
# import numpy as np

mol = Molecule("test/test_files/ethylene.xyz")


def test_Distance():
    """
    Test Distance reactivity functionality
    """
    d = reactivity.Distance(1, 2)
    s = d.get_settings(mol=mol)

    assert isinstance(s, Settings)


def test_Angle():
    """
    Test Angle reactivity functionality
    """
    ang = reactivity.Angle(1, 2, 3)
    s = ang.get_settings(mol=mol)

    assert isinstance(s, Settings)


def test_Dihedral():
    """
    Test Angle reactivity functionality
    """
    dihed = reactivity.Dihedral(1, 2, 3, 4)
    s = dihed.get_settings(mol=mol)

    assert isinstance(s, Settings)
