"""A module with named tuples (:class:`~typing.NamedTuple`)."""

from typing import Any, Type, NamedTuple, Callable, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from pyparsing import ParserElement
else:
    ParserElement = 'pyparsing.ParserElement'

__all__ = ['AtomBasisKey', 'AtomBasisData', 'AtomXYZ', 'CGF',
           'InfoMO', 'InputKey', 'MO', 'MO_metadata', 'ParseWarning']


class AtomBasisKey(NamedTuple):
    atom: Any
    basis: Any
    basisFormat: Any


class AtomBasisData(NamedTuple):
    exponents: Any
    coefficients: Any


class AtomXYZ(NamedTuple):
    symbol: Any
    xyz: Any


class CGF(NamedTuple):
    primitives: Any
    orbType: Any


class InfoMO(NamedTuple):
    eigenVals: Any
    coeffs: Any


class InputKey(NamedTuple):
    name: Any
    args: Any


class MO(NamedTuple):
    coordinates: Any
    cgfs: Any
    coefficients: Any


class MO_metadata(NamedTuple):
    """Molecular Orbitals Parsing.

    .. code::

        MO EIGENVALUES, MO OCCUPATION NUMBERS, AND SPHERICAL MO EIGENVECTORS

        1                      2
                                 -0.9857682741370732    -0.9831467097855797

                                  2.0000000000000000     2.0000000000000000

            1     1 cd  2s       -0.0015026981889089    -0.0103313715516893
            2     1 cd  3s       -0.0005376142747880    -0.0041729598190025
            3     1 cd  3py      -0.0013790317507575     0.0132729535025288
            4     1 cd  3pz      -0.0015557487597731    -0.0005486094359245
            5     1 cd  3px      -0.0013339995106232    -0.0100914249163043
            6     1 cd  4py      -0.0003884918433452     0.0046068283721132

    """

    nOccupied: Any
    nOrbitals: Any
    nOrbFuns: Any


class ParseWarning(NamedTuple):
    warn_type: Type[Warning]
    parser: ParserElement
    func: Callable[[str], Optional[str]] = NotImplemented

    @staticmethod
    def return_msg(msg: str) -> str:
        """Return the passed *msg* in unaltered form."""
        return msg


# Replace NotImplemented with ParseWarning.return_msg()
ParseWarning.__new__.__defaults__ = (ParseWarning.return_msg,)  # type: ignore
