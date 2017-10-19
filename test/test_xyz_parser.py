from nose.plugins.attrib import attr
from qmflows.parsers import (manyXYZ, parse_string_xyz)
from qmflows.utils import (chunksOf)


@attr('fast')
def test_multiple_geometries():
    """
    Test the reading of multiples molecular geometries from a file.
    """
    path_xyz = 'test/test_files/five_points_ethylene.xyz'

    with open(path_xyz, 'r') as f:
        ls = f.readlines()

    xs = [''.join(x) for x in chunksOf(ls, 8)]

    assert list(map(parse_string_xyz, xs)) == manyXYZ(path_xyz)
