__author__ = "Felipe Zapata"

__all__ = ["anyChar", "floatNumber", "floatNumberDot", "minusOrplus",
           "parse_file", "parse_section", "skipAnyChar"]

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
skipLine = Suppress(SkipTo('\n'))


# Generic Functions


def parse_file(p, file_name):
    """
    Wrapper over parseFile method
    """
    try:
        r = p.parseFile(file_name)
        return r[0]
    except ParseException:
        msg = "Error Trying to parse {}".format(p.name)
        print(msg)
        raise


def parse_section(start, end):
    """
    Read the lines from `start` to `end`.
    """
    s = Literal('{}'.format(start))
    e = Literal('{}'.format(end))

    return Suppress(SkipTo(s)) + skipLine + SkipTo(e)
