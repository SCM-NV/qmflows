from collections import namedtuple
from typing import Type, NamedTuple, Callable, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from pyparsing import ParserElement
else:
    ParserElement = 'pyparsing.ParserElement'

__all__ = ['AtomBasisKey', 'AtomBasisData', 'AtomXYZ', 'CGF',
           'InfoMO', 'InputKey', 'MO', 'ParseWarning']


def _return_msg(msg: str) -> str:
    return msg


# Named Tuples
AtomBasisKey = namedtuple("AtomBasisKey", ("atom", "basis", "basisFormat"))
AtomBasisData = namedtuple("AtomBasisData", ("exponents", "coefficients"))
AtomXYZ = namedtuple("AtomXYZ", ("symbol", "xyz"))
CGF = namedtuple("CGF", ("primitives", "orbType"))
InfoMO = namedtuple("InfoMO", ("eigenVals", "coeffs"))
InputKey = namedtuple("InpuKey", ("name", "args"))
MO = namedtuple("MO", ("coordinates", "cgfs", "coefficients"))


class ParseWarning(NamedTuple):
    warn_type: Type[Warning]
    parser: ParserElement
    func: Callable[[str], Optional[str]] = _return_msg
