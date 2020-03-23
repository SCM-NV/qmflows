"""Test the xyz file format readers."""
from assertionlib import assertion
from more_itertools import chunked

from qmflows.parsers import manyXYZ, parse_string_xyz, readXYZ
from qmflows.test_utils import PATH, PATH_MOLECULES


def test_multiple_geometries():
    """Test the reading of multiples molecular geometries from a file."""
    path_xyz = PATH / "molecules" / "five_points_ethylene.xyz"

    with open(path_xyz, 'r') as f:
        ls = f.readlines()

    xs = [''.join(x) for x in chunked(ls, 8)]

    assertion.eq(list(map(parse_string_xyz, xs)), manyXYZ(path_xyz))


def test_xyz_reader():
    """Test the xyz parser from a file."""
    path_qd = PATH_MOLECULES / "Cd68Cl26Se55__26_acetate.xyz"
    mol = readXYZ(path_qd)
    for atom in mol:
        assertion.le(len(atom.symbol), 2)
        assertion.len_eq(atom.xyz, 3)
