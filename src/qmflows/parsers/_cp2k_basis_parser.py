"""Functions for reading CP2K basis set files."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING
from itertools import islice

import numpy as np

from ..type_hints import PathLike
from ..common import AtomBasisData, AtomBasisKey

if TYPE_CHECKING or sys.version_info >= (3, 9):
    _Basis2Tuple = tuple[list[AtomBasisKey], list[AtomBasisData]]
    from collections.abc import Iterator, Iterable
else:
    from typing import Iterator

__all__ = ["read_cp2k_basis"]


class _BasisFileIter(Iterator[str]):
    """Enumerate through the passed ``iterable`` and remove all empty and commented lines."""

    __slots__ = ("__weakref__", "_enumerator", "_name", "_index")

    @property
    def index(self) -> int:
        """Get the index within the current iterator."""
        return self._index

    def __init__(self, iterable: Iterable[str], start: int = 0) -> None:
        self._enumerator = ((i, j.strip().rstrip()) for i, j in enumerate(iterable, start=start))
        self._name: "str | None" = getattr(iterable, "name", None)
        self._index = start

    def __iter__(self) -> _BasisFileIter:
        return self

    def __next__(self) -> str:
        while True:
            self._index, item = next(self._enumerator)
            if item and not item.startswith("#"):
                return item

    def __repr__(self) -> str:
        return f"<{type(self).__name__} name={self._name!r} index={self._index!r}>"


def _read_basis(f: _BasisFileIter, allow_multiple_exponents: bool) -> _Basis2Tuple:
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
            if n_sets != 1 and not allow_multiple_exponents:
                raise ValueError(
                    "Basis sets with more than 1 set of exponents "
                    "require `allow_multiple_exponens=True`"
                )

            for n in range(n_sets):
                # Parse the basis format, its exponents and its coefficients
                basis_fmt = tuple(int(j) for j in next(f).split())
                n_exp = basis_fmt[3]
                basis_data = np.array([j.split() for j in islice(f, 0, n_exp)], dtype=np.float64)
                exp, coef = basis_data[:, 0], basis_data[:, 1:].T
                for basis in basis_list:
                    keys.append(AtomBasisKey(atom, basis, basis_fmt, n))
                    values.append(AtomBasisData(exp, coef))
        except Exception as ex:
            basis = f"{atom} {basis_list[0]}" if len(basis_list) > 0 else atom
            raise ValueError(f'Failed to parse the {basis!r} basis set on line {f.index}') from ex
    return keys, values


def read_cp2k_basis(file: PathLike, *, allow_multiple_exponents: bool = False) -> _Basis2Tuple:
    """Read the Contracted Gauss function primitives format from a text file.

    Parameters
    ----------
    file : path-like
        The to-be read CP2K basis set file.
    allow_multiple_exponents : bool
        Whether to allow the parsing of basis sets consisting of multiple exponent sets.

    """
    with open(file, "r") as f:
        return _read_basis(_BasisFileIter(f, start=1), allow_multiple_exponents)
