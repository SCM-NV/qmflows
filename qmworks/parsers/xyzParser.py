__author__ = "Felipe Zapata"

__all__ = ['parse_string_xyz', 'readXYZ', 'manyXYZ']

# ==========> Standard libraries and third-party <===============
from collections import namedtuple
from pyparsing   import *

# ==================> Internal modules <====================
from qmworks.parsers.parser import floatNumber, natural
from qmworks.utils import zipWith
# ================================================================================
AtomXYZ = namedtuple("AtomXYZ", ("symbol", "xyz"))

header  = natural + LineEnd() + restOfLine

label   = Word(alphas, max=2)

xyz     = floatNumber * 3

atomParser = label.setResultsName("label") + xyz.setResultsName("xyz")

parser_xyz = Suppress(header) + OneOrMore(Group(atomParser))

# ==============================<>====================================


def parse_string_xyz(xs):
    """
    :param: xs
    :type:  string
    """
    rs = parser_xyz.parseString(xs)
    return createAtoms(rs)


def readXYZ(pathXYZ):
    """
    Parse From File
    :param: pathXYZ
    :type:  string
    """
    xs = parser_xyz.parseFile(pathXYZ)
    return createAtoms(xs)


def manyXYZ(pathXYZ):
    """
    :param: pathXYZ
    :type:  string
    """
    manyMol = OneOrMore(Group(parser_xyz))
    xss     = manyMol.parseFile(pathXYZ)
    return list(map(createAtoms, xss))


def createAtoms(xs):
    ls = [a.label.lower() for a in xs]
    rs = [list(map(float, a.xyz)) for a in xs]
    return zipWith(AtomXYZ)(ls)(rs)
