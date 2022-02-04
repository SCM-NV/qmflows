"""Functions for reading CP2K basis set files."""

from typing import Iterator, Generator, Tuple, List
from itertools import islice

import numpy as np

from ..type_hints import PathLike
from ..common import AtomBasisData, AtomBasisKey

__all__ = ["readCp2KBasis"]

_Basis2Tuple = Tuple[List[AtomBasisKey], List[AtomBasisData]]


def _basis_file_iter(f: Iterator[str]) -> Generator[str, None, None]:
    """Iterate through `f` and remove all empty and commented lines."""
    for i in f:
        i = i.strip().rstrip("\n")
        if not i or i.startswith("#"):
            continue
        yield i


def _read_basis(iterator: Iterator[str]) -> _Basis2Tuple:
    """Helper function for parsing the opened basis set file."""
    f = _basis_file_iter(iterator)
    keys = []
    values = []
    for i in f:
        # Read the atom type and basis set name(s)
        atom, *basis_list = i.split()
        atom = atom.capitalize()

        # Identify the number of exponent sets
        n_sets = int(next(f))
        if n_sets != 1:
            raise NotImplementedError(
                "Basis sets with more than 1 set of exponents are not supported yet"
            )

        for _ in range(n_sets):
            # Parse the basis format, its exponents and its coefficients
            basis_fmt = [int(j) for j in next(f).split()]
            n_exp = basis_fmt[3]
            basis_data = np.array([j.split() for j in islice(f, 0, n_exp)], dtype=np.float64)
            exp, coef = basis_data[:, 0], basis_data[:, 1:]

            # Two things happen whenever an basis set alias is encountered (i.e. `is_alias > 0`):
            # 1. The `alias` field is set for the keys
            # 2. The `AtomBasisData` instance, used for the original value, is reused
            for is_alias, basis in enumerate(basis_list):
                if not is_alias:
                    basis_key = AtomBasisKey(atom, basis, basis_fmt)
                    basis_value = AtomBasisData(exp, coef)
                    keys.append(basis_key)
                else:
                    keys.append(AtomBasisKey(atom, basis, basis_fmt, alias=basis_key))
                values.append(basis_value)
    return keys, values


def readCp2KBasis(file: PathLike) -> _Basis2Tuple:
    """Read the Contracted Gauss function primitives format from a text file."""
    with open(file, "r") as f:
        return _read_basis(f)
