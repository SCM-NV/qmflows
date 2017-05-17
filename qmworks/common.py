
__author__ = "Felipe Zapata"

__all__ = ['AtomBasisKey', 'AtomBasisData', 'AtomXYZ', 'CGF',
           'InfoMO', 'InputKey', 'MO']


# ================> Python Standard  and third-party <==========
from collections import namedtuple

# ======================================================================
# Named Tuples
AtomBasisKey = namedtuple("AtomBasisKey", ("atom", "basis", "basisFormat"))
AtomBasisData = namedtuple("AtomBasisData", ("exponents", "coefficients"))
AtomXYZ = namedtuple("AtomXYZ", ("symbol", "xyz"))
CGF = namedtuple("CGF", ("primitives", "orbType"))
InfoMO = namedtuple("InfoMO", ("eigenVals", "coeffs"))
InputKey = namedtuple("InpuKey", ("name", "args"))
MO = namedtuple("MO", ("coordinates", "cgfs", "coefficients"))
