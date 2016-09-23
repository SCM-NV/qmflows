
__all__ = ['parse_dipole', 'parse_frequencies', 'parse_hessian',
           'parse_gradient']


from pyparsing import (alphanums, Literal, OneOrMore, SkipTo, Suppress,
                       Word)

from .parser import (floatNumber, integer, parse_file, parse_section,
                     skipSupress)

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


def parse_dipole(file_name):
    """
    Parse dipole moment from the *.dat file.
    """
    l = Literal('DIPOLE')
    p = Suppress(SkipTo(l) + l) + OneOrMore(floatNumber)
    rs = parse_file(p, file_name).asList()
    return string_array_to_float(rs)


def parse_gradient(file_name):
    """
    Parse Gradient from the *.dat file.
    """
    p = Suppress(Word(alphanums) + floatNumber) + OneOrMore(floatNumber)
    return parse_dollar_section(file_name, p, '$GRAD')


def parse_hessian(file_name):
    """
    Parse the hessian from the *.dat produced by gamess.
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

# def parse_energy():
#     """
#     parse Total energy.
#     """
#     Suppress(SkipTo(Regex(r'^E'))
