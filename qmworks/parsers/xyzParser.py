__author__ = "Felipe Zapata"

__all__ = ['parse_string_xyz', 'readXYZ', 'manyXYZ']

# ===================> Standard libraries and third-party <====================
from pyparsing   import (alphas, Group, LineEnd, OneOrMore, restOfLine,
                         Suppress, Word)

from qmworks.common import AtomXYZ
from qmworks.parsers.parser import (floatNumber, natural)
from qmworks.utils import zipWith
# =============================================================================

header  = natural + LineEnd() + restOfLine

label   = Word(alphas, max=2)

xyz     = floatNumber * 3

atomParser = label.setResultsName("label") + xyz.setResultsName("xyz")

parser_xyz = Suppress(header) + OneOrMore(Group(atomParser))

# ==============================<>====================================


def parse_string_xyz(xs):
    """
    Read a molecula geometry in XYZ format from a string.

    :param: xs
    :type:  string
    :return: [AtomXYZ]
    """
    rs = parser_xyz.parseString(xs)
    return createAtoms(rs)


def readXYZ(pathXYZ):
    """
    Parse molecular geometry in XYZ format from a file.

    :param: pathXYZ
    :type:  string
    :return: [AtomXYZ]
    """
    xs = parser_xyz.parseFile(pathXYZ)
    return createAtoms(xs)


def manyXYZ(pathXYZ):
    """
    Read one or more molecular geometries in XYZ format from a file.

    :param: pathXYZ
    :type:  string
    :return: [[AtomXYZ]]
    """
    manyMol = OneOrMore(Group(parser_xyz))
    xss     = manyMol.parseFile(pathXYZ)
    return list(map(createAtoms, xss))


def createAtoms(xs):
    """
    Create an AtomXYZ tuple from a string.
    """
    ls = [a.label.lower() for a in xs]
    rs = [list(map(float, a.xyz)) for a in xs]
    return zipWith(AtomXYZ)(ls)(rs)
