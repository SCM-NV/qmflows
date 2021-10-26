"""common (:class:`~typing.NamedTuple`) used mostly for parsing."""
from typing import (Any, Type, NamedTuple, Callable, Optional, Tuple,
                    Sequence, Union, TYPE_CHECKING)

import numpy as np

if TYPE_CHECKING:
    from pyparsing import ParserElement
else:
    ParserElement = 'pyparsing.ParserElement'

__all__ = ['AtomBasisKey', 'AtomBasisData', 'AtomXYZ', 'CGF',
           'InfoMO', 'MO_metadata', 'ParseWarning', 'CP2KVersion']


class AtomBasisKey(NamedTuple):
    """Namedtuple containing the `basisFormat` for a given `basis` and `atom`."""

    atom: str
    basis: str
    basisFormat: Union[Sequence[int],
                       Sequence[Tuple[str, int]]]
    # Orca uses 2-tuples while CP2K uses integer


class AtomBasisData(NamedTuple):
    """Contraction coefficients and exponents for a gaussian basis."""

    exponents: np.ndarray
    coefficients: np.ndarray


class AtomXYZ(NamedTuple):
    """Symbol and coordinates of an Atom."""

    symbol: str
    xyz: Tuple[float, float, float]


class CGF(NamedTuple):
    """Contracted Gaussian function."""

    primitives: tuple  # (AtomBasisKey, AtomBasisData)
    orbType: str  # (S, Px, Py, Pz, etc.)


class InfoMO(NamedTuple):
    """Energies and coefficients of the molecular orbitals."""

    eigenvalues: np.ndarray  # Orbitals eigenvalues
    eigenvectors: np.ndarray  # Orbitals eigenvectors


class MO_metadata(NamedTuple):
    """Metadata for the molecular orbitals.

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

    nOccupied: int
    nOrbitals: int
    nOrbFuns: int
    nspinstates: int = 1


class ParseWarning(NamedTuple):
    """Generate a Namedtuple using a sort of data class.

    see: https://docs.python.org/3/library/typing.html#typing.NamedTuple
    """

    warn_type: Type[Warning]
    parser: ParserElement
    func: Callable[[str], Optional[str]] = NotImplemented

    @staticmethod
    def return_msg(msg: str) -> str:
        """Return the passed *msg* in unaltered form."""
        return msg


# Replace NotImplemented with ParseWarning.return_msg()
ParseWarning.__new__.__defaults__ = (ParseWarning.return_msg,)  # type: ignore


class CP2KVersion(NamedTuple):
    """A named 2-tuple with the CP2K major and minor version."""

    major: int
    minor: int
