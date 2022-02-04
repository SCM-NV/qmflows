"""Functions for reading CP2K basis set files."""

from typing import Iterator, Tuple, List, Tuple, Iterable
from itertools import islice

import numpy as np

from ..type_hints import PathLike
from ..common import AtomBasisData, AtomBasisKey

__all__ = ["readCp2KBasis"]

_Basis2Tuple = Tuple[List[AtomBasisKey], List[AtomBasisData]]


class _BasisFileIter:
    """Enumerate through the passed ``iterable`` and remove all empty and commented lines."""

    __slots__ = ("__weakref__", "_enumerator", "_name", "_index")

    @property
    def index(self) -> int:
        """Get the index within the current iterator."""
        return self._index

    def __init__(self, iterable: Iterable[str], start: int = 0) -> None:
        self._enumerator = enumerate(iterable, start=start)
        self._name: "str | None" = getattr(iterable, "name", None)
        self._index = start

    def __iter__(self) -> "_BasisFileIter":
        return self

    def __next__(self) -> str:
        i, item = next(self._enumerator)
        item = item.strip().rstrip("\n")
        if not item or item.startswith("#"):
            return self.__next__()

        self._index = i
        return item

    def __repr__(self) -> str:
        return f"<{type(self).__name__} name={self._name!r} index={self._index!r}>"


def _read_basis(f: _BasisFileIter) -> _Basis2Tuple:
    """Helper function for parsing the opened basis set file."""
    keys = []
    values = []
    for i in f:
        # Read the atom type and basis set name(s)
        atom, *basis_list = i.split()
        atom = atom.capitalize()

        try:
            # Identify the number of exponent sets
            n_sets = int(next(f))
            if n_sets != 1:
                raise NotImplementedError(
                    "Basis sets with more than 1 set of exponents are not supported yet"
                )

            for _ in range(n_sets):
                # Parse the basis format, its exponents and its coefficients
                basis_fmt = tuple(int(j) for j in next(f).split())
                n_exp = basis_fmt[3]
                basis_data = np.array([j.split() for j in islice(f, 0, n_exp)], dtype=np.float64)
                exp, coef = basis_data[:, 0], basis_data[:, 1:].T
                for basis in basis_list:
                    keys.append(AtomBasisKey(atom, basis, basis_fmt))
                    values.append(AtomBasisData(exp, coef))
        except Exception as ex:
            raise ValueError(
                f'Failed to parse the basis set "{atom} {basis_list[0]}" on line {f.index}'
            ) from ex
    return keys, values


def readCp2KBasis(file: PathLike) -> _Basis2Tuple:
    """Read the Contracted Gauss function primitives format from a text file."""
    with open(file, "r") as f:
        return _read_basis(_BasisFileIter(f, start=1))
