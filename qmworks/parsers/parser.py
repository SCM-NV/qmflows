__author__ = "Felipe Zapata"

# ===============> Standard libraries and third-party <========================
from pyparsing import (CaselessKeyword, Combine, Literal, nums, Optional,
                       ParseException, Regex, SkipTo, Suppress, Word)

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
