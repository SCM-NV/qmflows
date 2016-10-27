__author__ = "Felipe Zapata"

__all__ = ['parse_string_xyz', 'readXYZ', 'manyXYZ', "string_to_plams_Molecule"]

# ===================> Standard libraries and third-party <====================
from plams import (Atom, Molecule)
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


def tuplesXYZ_to_plams(xs):
    """ Transform a list of namedTuples to a Plams molecule """
    plams_mol = Molecule()
    for at in xs:
        symb = at.symbol
        cs = at.xyz
        plams_mol.add_atom(Atom(symbol=symb, coords=tuple(cs)))

    return plams_mol


def string_to_plams_Molecule(xs):
    """Convert a molecule stored in a string to a plams Molecule"""
    return  tuplesXYZ_to_plams(parse_string_xyz(xs))
