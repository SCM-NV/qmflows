"""common (:class:`~typing.NamedTuple`) used mostly for parsing."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any, NamedTuple, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from numpy import float64 as f8, int64 as i8
    import pyparsing

__all__ = ['AtomBasisKey', 'AtomBasisData', 'AtomXYZ', 'CGF', 'CP2KInfoMO',
           'InfoMO', 'MO_metadata', 'ParseWarning', 'CP2KVersion']


class AtomBasisKey(NamedTuple):
    """Namedtuple containing the `basisFormat` for a given `basis` and `atom`."""

    #: Atomic symbol.
    atom: str

    #: Name of the basis set.
    basis: str

    # NOTE: Orca uses 2-tuples while CP2K uses integer
    #: The basis set format.
    basisFormat: Sequence[int] | Sequence[tuple[str, int]]

    #: The index of the exponent set.
    #: Relevant for basis sets consisting of multiple :attr:`basisFormat` blocks.
    exponent_set: int = 0


class AtomBasisData(NamedTuple):
    """Contraction coefficients and exponents for a gaussian basis."""

    exponents: NDArray[f8]
    coefficients: NDArray[f8]


class AtomXYZ(NamedTuple):
    """Symbol and coordinates of an Atom."""

    symbol: str
    xyz: tuple[float, float, float]


class CGF(NamedTuple):
    """Contracted Gaussian function."""

    primitives: tuple  # (AtomBasisKey, AtomBasisData)
    orbType: str  # (S, Px, Py, Pz, etc.)


class InfoMO(NamedTuple):
    """Energies and coefficients of the molecular orbitals."""

    #: Orbital eigenvalues.
    eigenvalues: NDArray[f8]

    #: Orbital eigenvectors.
    eigenvectors: NDArray[f8]


class CP2KInfoMO(NamedTuple):
    """Energies and coefficients of the CP2K molecular orbitals."""

    #: Orbital eigenvalues.
    eigenvalues: NDArray[f8]

    #: Orbital eigenvectors.
    eigenvectors: NDArray[f8]

    #: MO indices.
    orb_index: NDArray[i8]

    #: Occupation numbers.
    occupation: NDArray[f8]

    def get_nocc_nvirt(self, threshold: None | float = None) -> tuple[int, int]:
        """Return the number of occupied and virtual orbitals within the MO range spanned \
        by this instance.

        Parameters
        ----------
        threshold : None | float
            The occupation number threshold for determining what consitutes an occupied orbital.
            If ``None`` assume that all occupied orbitals are defined by a non-zero occupation
            number (*e.g.* canonical orbitals).

        Returns
        -------
        tuple[int, int]
            A 2-tuple with the number of occupied and virtual orbitals.

        """
        if threshold is None:
            is_occ = self.occupation != 0
        else:
            is_occ = self.occupation >= threshold
        nocc = is_occ.sum().item()
        return nocc, len(self.occupation) - nocc

    def __array__(self, dtype: None | np.dtype[Any] = None) -> NDArray[Any]:
        """Convert this instance into a structured array."""
        struc_dtype = np.dtype([
            ("eigenvalues", "f8"),
            ("eigenvectors", "f8", (self.eigenvectors.shape[0],)),
            ("orb_index", "i8"),
            ("occupation", "f8"),
        ])

        ret = np.empty(self.eigenvalues.shape[0], dtype=struc_dtype)
        ret["eigenvalues"] = self.eigenvalues
        ret["eigenvectors"] = np.abs(self.eigenvectors.T)
        ret["orb_index"] = self.orb_index
        ret["occupation"] = self.occupation
        return ret if dtype is None else ret.astype(dtype)


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

    #: The number of occupied orbitals.
    #: Has either 1 or 2 elements depending on whether they're spin-orbitals or not.
    nOccupied: list[int]

    #: The number of orbitals.
    #: Has either 1 or 2 elements depending on whether they're spin-orbitals or not.
    nOrbitals: list[int]

    #: The number of basis functions.
    nOrbFuns: int

    @property
    def nspinstates(self) -> int:
        """The number of spin states."""
        return len(self.nOrbitals)


class ParseWarning(NamedTuple):
    """Generate a Namedtuple using a sort of data class.

    see: https://docs.python.org/3/library/typing.html#typing.NamedTuple
    """

    warn_type: type[Warning]
    parser: pyparsing.ParserElement
    func: Callable[[str], None | str] = NotImplemented

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
