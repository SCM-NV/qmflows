__author__ = "Felipe Zapata"

# ===============> Standard libraries and third-party <========================
from plams import (Atom, Molecule)
from pyparsing import (CaselessKeyword, Combine, Literal, nums, Optional,
                       ParseException, Regex, SkipTo, Suppress, Word)
import numpy as np

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
anyChar     = Regex('.')
skipAnyChar = Suppress(anyChar)
skipSupress = lambda z: Suppress(SkipTo(z))
skipLine = Suppress(skipSupress('\n'))


# Generic Functions

def parse_file(p, file_name):
    """
    Wrapper over the parseFile method
    """
    try:
        return p.parseFile(file_name)
    except ParseException:
        msg = "Error Trying to parse: {} in file: {}".format(p, file_name)
        print(msg)
        raise


def parse_section(start, end):
    """
    Read the lines from `start` to `end`.
    """
    s = Literal('{}'.format(start))
    e = Literal('{}'.format(end))

    return Suppress(SkipTo(s)) + skipLine + SkipTo(e)


def string_array_to_molecule(parser_fun, file_name, mol=None):
    """
    Convert a Numpy string array like:

    [['C', '-1.487460', '-0.028670', '-0.000060'],
    ['O', '0.376340', '0.028670', '-0.000060'],
    ['H', '-1.818910', '-1.067060', '-0.000060'],
    ['H', '-1.866470', '0.473700', '0.889930'],
    ['H', '-1.866470', '0.473700', '-0.890040'],
    ['H', '0.756720', '-0.950010', '-0.000060']]

    To a plams ``Molecule``.
    """
    string_array_to_float = np.vectorize(float)
    mols = parse_file(parser_fun, file_name).asList()
    last_mol = np.array(mols[-1])
    elems = last_mol[:, 0]
    coords = string_array_to_float(last_mol[:, 1:])
    if mol:
        if len(coords) == len(mol):
            plams_mol = mol
            for i in range(len(plams_mol)):
                plams_mol.atoms[i].coords = tuple([float(c) for c in coords[i]])
        else:
            raise RuntimeError('Output molecule does not match input molecule')
    else:
        plams_mol = Molecule()
        for e, c in zip(elems, coords):
            plams_mol.add_atom(Atom(symbol=e, coords=tuple(c)))
    return plams_mol
