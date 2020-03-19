"""General utilities to parse both input/out files."""

__all__ = ['anyChar', 'integer', 'natural', 'parse_file', 'parse_section',
           'skipAnyChar', 'skipLine', 'skipSupress', 'try_search_pattern']


import re
import typing

import numpy as np
from ..type_hints import PathLike
from pyparsing import (CaselessKeyword, Combine, Literal, Optional,
                       ParseException, ParserElement, ParseResults, Regex,
                       SkipTo, Suppress, Word, nums)
from scm.plams import Atom, Molecule

# Literals
point = Literal('.')
e = CaselessKeyword('E')
minusOrplus = Literal('+') | Literal('-')

# Parsing Floats
natural = Word(nums)
integer = Combine(Optional(minusOrplus) + natural)
floatNumber = Regex(r'(\-)?\d+(\.)(\d*)?([eE][\-\+]\d+)?')

floatNumberDot = Regex(r'(\-)?(\d+)?(\.)(\d*)?([eE][\-\+]\d+)?')

# Parse Utilities


def skipSupress(z: str) -> ParserElement:
    """Skip until `z` and suppress the skipped values."""
    return Suppress(SkipTo(z))


anyChar = Regex('.')
skipAnyChar = Suppress(anyChar)
skipLine = Suppress(skipSupress('\n'))

# Generic Functions


def parse_file(p: ParserElement, file_name: PathLike) -> ParseResults:
    """Apply parser `p` on file `file_name`."""
    try:
        return p.parseFile(file_name)
    except ParseException:
        msg = "Error Trying to parse: {} in file: {}".format(p, file_name)
        print(msg)
        raise


def parse_section(start: str, end: str) -> ParserElement:
    """Read the lines from `start` to `end`."""
    s = Literal('{}'.format(start))
    e = Literal('{}'.format(end))

    return Suppress(SkipTo(s)) + skipLine + SkipTo(e)


def string_array_to_molecule(
        parser_fun: ParserElement, file_name: PathLike, mol: Molecule = None) -> Molecule:
    """Convert a Numpy string array.

    It takes an array like:
    [['C', '-1.487460', '-0.028670', '-0.000060'],
    ['O', '0.376340', '0.028670', '-0.000060'],
    ['H', '-1.818910', '-1.067060', '-0.000060'],
    ['H', '-1.866470', '0.473700', '0.889930'],
    ['H', '-1.866470', '0.473700', '-0.890040'],
    ['H', '0.756720', '-0.950010', '-0.000060']]

    and covert it to a plams ``Molecule``.
    """
    mols = parse_file(parser_fun, file_name).asList()
    last_mol = np.array(mols[-1])
    elems = last_mol[:, 0]
    coords = np.array(last_mol[:, 1:], dtype=float)
    if mol:
        if len(coords) == len(mol):
            mol.from_array(coords)
        else:
            raise RuntimeError('Output molecule does not match input molecule')
    else:
        plams_mol = Molecule()
        for e, c in zip(elems, coords):
            plams_mol.add_atom(Atom(symbol=e, coords=tuple(c)))
    return plams_mol


def try_search_pattern(pat: str, file_name: PathLike) -> typing.Optional[str]:
    """Search for an specific pattern in  a file."""
    try:
        with open(file_name, 'r') as f:
            for line in f:
                if re.search(pat, line):
                    return line
    except NameError:
        return None
    except FileNotFoundError:
        msg2 = 'There is not a file: {}\n'.format(file_name)
        raise RuntimeError(msg2)
