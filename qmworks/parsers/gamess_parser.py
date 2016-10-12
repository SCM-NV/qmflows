
__all__ = ['parse_dipole', 'parse_total_energy', 'parse_frequencies',
           'parse_gradient', 'parse_hessian', 'parse_molecule']

from pyparsing import (alphanums, Group, Literal, OneOrMore, Regex, Suppress,
                       Word)
from .parser import (floatNumber, integer, parse_file, parse_section, skipLine,
                     skipSupress, string_array_to_molecule)

import numpy as np

string_array_to_float = np.vectorize(float)


def parse_dollar_section(file_name, parser, start, trailing_lines=1):
    """
    Parse a Section starting with a $Header and ending with the $END token.
    """
    xs = parse_file(parse_section(start, '$END'), file_name)[0]
    # Skip first and last lines
    lines = xs.splitlines()[trailing_lines:-1]
    xss = [parser.parseString(x).asList() for x in lines]

    return string_array_to_float(np.concatenate(xss))


def parse_molecule(file_name):
    """
    Parse the last available geometry from both the job_name.dat or job_name.log.
    """
    # Header in .out
    header1 = Literal("COORDINATES OF ALL ATOMS ARE (ANGS)")
    # Header in .dat
    header2 = Literal("COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)")
    p1 = skipSupress(header1 | header2) + \
        skipLine
    p2 = skipSupress("---\n") + skipLine
    parse_atoms = Group(Word(alphanums) + Suppress(floatNumber) +
                        floatNumber * 3)
    parse_mol = p1 + p2 + Group(OneOrMore(parse_atoms))
    many_mols = OneOrMore(parse_mol)

    return string_array_to_molecule(many_mols, file_name)


def parse_dipole(file_name):
    """
    Parse dipole moment from the output.dat file.
    """
    l = Literal('DIPOLE')
    p = skipSupress(l) + OneOrMore(floatNumber)
    rs = parse_file(p, file_name).asList()
    return string_array_to_float(rs)


def parse_gradient(file_name):
    """
    Parse Gradient from the job_name.dat file.
    """
    p = Suppress(Word(alphanums) + floatNumber) + OneOrMore(floatNumber)
    return parse_dollar_section(file_name, p, '$GRAD')


def parse_hessian(file_name):
    """
    Parse the hessian from the job_name.dat produced by gamess.
    """
    p = Suppress(integer * 2) + OneOrMore(floatNumber)
    return parse_dollar_section(file_name, p, '$HESS')


def parse_frequencies(file_name):
    """
    Parse frequencies and normal modes
    """
    start = "START OF NORMAL MODES"
    end = "END OF NORMAL MODES"
    xs = parse_file(parse_section(start, end), file_name)[0]
    # Split the string in lines containing both the frequencies and the modes
    modes = xs.split('FREQUENCY')[1:]
    freq_parser = Suppress(Literal('=')) + OneOrMore(floatNumber)
    modes_parser = skipSupress(floatNumber) + OneOrMore(floatNumber)
    reader = freq_parser + modes_parser

    rs = [reader.parseString(x).asList() for x in modes]
    arr = string_array_to_float(np.array(rs))
    # Frequencies are the first element of each row
    frequencies = arr[:, 0]
    array_modes = arr[:, 1:]
    return frequencies, array_modes


def parse_total_energy(file_name):
    """
    Parse Energy from the output file.
    """
    skip = Literal('TOTAL ENERGY') + Regex('=')
    p = skipSupress(skip) + Suppress(skip) + floatNumber

    return float(parse_file(p, file_name).asList()[0])
