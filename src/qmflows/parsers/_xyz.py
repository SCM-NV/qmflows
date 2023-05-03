"""XYZ file format readers."""

from __future__ import annotations

import os
from collections.abc import Iterable

from pyparsing import (Group, LineEnd, OneOrMore, Suppress, Word, alphas,
                       restOfLine, ParseResults)
from scm.plams import Atom, Molecule

from .utils import floatNumber, natural
from ..common import AtomXYZ

__all__ = ['parse_string_xyz', 'readXYZ',
           'manyXYZ', "string_to_plams_Molecule"]

# =============================================================================

header = natural + LineEnd() + restOfLine

label = Word(alphas, max=2)

xyz = floatNumber * 3

atomParser = label.setResultsName("label") + xyz.setResultsName("xyz")

parser_xyz = Suppress(header) + OneOrMore(Group(atomParser))

# ==============================<>====================================


def parse_string_xyz(xs: str) -> list[AtomXYZ]:
    """Read a molecula geometry in XYZ format from a string.

    :param: xs
    :type:  string
    :return: [AtomXYZ]
    """
    rs = parser_xyz.parseString(xs)
    return createAtoms(rs)


def readXYZ(pathXYZ: str | os.PathLike[str]) -> list[AtomXYZ]:
    """Parse molecular geometry in XYZ format from a file.

    :param: pathXYZ
    :type:  string
    :return: [AtomXYZ]
    """
    xs = parser_xyz.parseFile(os.fspath(pathXYZ))
    return createAtoms(xs)


def manyXYZ(pathXYZ: str | os.PathLike[str]) -> list[list[AtomXYZ]]:
    """Read one or more molecular geometries in XYZ format from a file.

    :param: pathXYZ
    :type:  string
    :return: [[AtomXYZ]]
    """
    manyMol = OneOrMore(Group(parser_xyz))
    xss = manyMol.parseFile(os.fspath(pathXYZ))
    return [createAtoms(i) for i in xss]


def createAtoms(xs: ParseResults) -> list[AtomXYZ]:
    """Create an AtomXYZ tuple from a string."""
    ls = (a.label.lower() for a in xs)
    rs = (tuple(map(float, a.xyz)) for a in xs)
    return [AtomXYZ(symbol, xyz) for symbol, xyz in zip(ls, rs)]  # type: ignore


def tuplesXYZ_to_plams(xs: Iterable[AtomXYZ]) -> Molecule:
    """Transform a list of namedTuples to a Plams molecule."""
    plams_mol = Molecule()
    for at in xs:
        symb = at.symbol
        cs = at.xyz
        plams_mol.add_atom(Atom(symbol=symb, coords=tuple(cs)))

    return plams_mol


def string_to_plams_Molecule(xs: str) -> Molecule:
    """Convert a molecule stored in a string to a plams Molecule."""
    return tuplesXYZ_to_plams(parse_string_xyz(xs))
