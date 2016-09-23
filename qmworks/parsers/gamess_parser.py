
from pyparsing import (OneOrMore, Suppress)

from .parser import (floatNumber, integer, parse_file, parse_section)

import numpy as np


def parse_hess(file_name):
    """
    Parse the hessian from the *.dat produced by gamess.
    """
    # Skipt to Integer index and read the floats

    p = Suppress(integer * 2) + OneOrMore(floatNumber)
    xs = parse_file(parse_section('$HESS', '$END'), file_name)
    # Skip first and last lines
    lines = xs.splitlines()[1:-1]
    xss = [p.parseString(x).asList() for x in lines]
    hess = np.concatenate(list(xss))
    floatArray = np.vectorize(float)

    return floatArray(hess)


def parse_frequencies(file_name):
    """
    Parse frequencies and normal modes
    """
    start = "START OF NORMAL MODES"
    end = "END OF NORMAL MODES"
    xs = parse_file(parse_section(start, end), file_name)

# def parse_energy():
#     """
#     parse Total energy.
#     """
#     Suppress(SkipTo(Regex(r'^E'))
