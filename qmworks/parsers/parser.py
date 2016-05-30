__author__ = "Felipe Zapata"

__all__ = ["anyChar", "floatNumber", "floatNumberDot", "minusOrplus", "skipAnyChar"]

# ==========> Standard libraries and third-party <===============
from pyparsing import (CaselessKeyword, Combine, Literal,
                       nums, Optional, Regex, Suppress, Word)

# ======================================================================

# ========= Literals ==========
point = Literal('.')
e = CaselessKeyword('E')
minusOrplus = Literal('+') | Literal('-')

# ======= Parsing Floats =================
natural = Word(nums)
integer = Combine(Optional(minusOrplus) + natural)
# floatNumber = Combine( integer + point  +
#                        Optional(number) +
#                        Optional( e + integer )
#                      )

# floatNumber = Regex(r'(\-)?\d+(\.\d*)?([eE]\d+)?')

floatNumber = Regex(r'(\-)?\d+(\.)(\d*)?([eE][\-\+]\d+)?')

floatNumberDot = Regex(r'(\-)?(\d+)?(\.)(\d*)?([eE][\-\+]\d+)?')


# Parse Utilities

anyChar     = Regex('.')

skipAnyChar = Suppress(anyChar)


