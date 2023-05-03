"""Read Orca output files."""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
from more_itertools import chunked
from pyparsing import Group, OneOrMore, Word, alphanums
from scm.plams import Atom, Molecule

from .utils import (floatNumber, parse_file, parse_section, skipLine,
                    skipSupress, string_array_to_molecule,
                    try_search_pattern)
from ._xyz import manyXYZ
from ..common import InfoMO

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from numpy import float64 as f8

__all__ = [
    'parse_hessian', 'parse_frequencies', 'parse_molecule',
    'parse_molecular_orbitals', 'parse_molecule_traj', 'parse_normal_modes']


def parse_molecule(file_name: str | os.PathLike[str], mol: None | Molecule = None) -> Molecule:
    """Parse The Cartesian coordinates from the output file."""
    header = "CARTESIAN COORDINATES (ANGSTROEM)"
    p1 = skipSupress(header) + skipLine * 2
    parse_atoms = Group(Word(alphanums) + floatNumber * 3)
    parse_mol = p1 + Group(OneOrMore(parse_atoms))

    parse_many_mol = OneOrMore(parse_mol)

    return string_array_to_molecule(parse_many_mol, file_name, mol=mol)


def parse_molecule_traj(file_traj: str | os.PathLike[str]) -> Molecule:
    """Read Molecules from the job_name.traj file."""
    mols = manyXYZ(file_traj)
    # Last geometry corresponds to the optimized structure
    opt_mol = mols[-1]

    plams_mol = Molecule()
    for at in opt_mol:
        symb = at.symbol
        cs = at.xyz
        plams_mol.add_atom(Atom(symbol=symb, coords=tuple(cs)))

    return plams_mol


def parse_hessian(file_hess: str | os.PathLike[str], start: str = '$hessian') -> NDArray[f8]:
    """Read the hessian matrix in cartesian coordinates from the job_name.hess file.

    :returns: Numpy array
    """
    return read_blocks_from_file(start, '\n\n', file_hess)


def parse_normal_modes(file_hess: str | os.PathLike[str]) -> NDArray[f8]:
    """Return the normal modes from the job_name.hess file."""
    start = '$normal_modes'
    return read_blocks_from_file(start, '\n\n', file_hess)


def parse_frequencies(file_hess: str | os.PathLike[str]) -> NDArray[f8]:
    """Parse the vibrational frequencies from the job_name.hess file."""
    p = parse_section('$vibrational_frequencies', '\n\n')
    lines = parse_file(p, file_hess)[0].splitlines()

    return np.array([x.split()[-1] for x in lines[1:]], dtype=np.float64)


def read_blocks_from_file(start: str, end: str, file_name: str | os.PathLike[str]) -> NDArray[f8]:
    """Read a matrix printed in block format.

    :param  start: token identifying the start of the block.
    :param end: characters signaling end of block.
    :param file_name: Name of the file containing the matrix printed in blocks
    :returns: Numpy array
    """
    p = parse_section(start, end)
    raw = parse_file(p, file_name)[0].splitlines()
    number_of_basis = int(raw[0].split()[0])
    lines = raw[1:]

    # Matrix elements are printed in block of 6 columns
    nblocks = number_of_basis // 6
    rest = number_of_basis % 6
    nblocks = nblocks if rest == 0 else nblocks + 1

    # a block start with a line header then the data
    blocks = [read_block(block)
              for block in chunked(lines, number_of_basis + 1)]

    return np.concatenate(blocks, axis=1)


def read_block(lines: Sequence[str]) -> NDArray[f8]:
    """Read a block containing the values of the block matrix.

    The format is similar to:

              0          1          2          3          4          5
      0   0.752994  0.078644   0.000000  -0.134007   0.156950  -0.000000
      1   0.078857  0.848881  -0.000000   0.162989  -0.312155   0.000000
      2   0.000000 -0.000000   0.074907  -0.000000   0.000000  -0.026031
      3  -0.133981  0.162928  -0.000000   0.133434  -0.164602   0.000000

    Read the matrix skiping the header and the first integer index

    :returns: Numpy array
    """
    return np.array([x.split()[1:] for x in lines[1:]], dtype=np.float64)


def parse_molecular_orbitals(file_name: str | os.PathLike[str]) -> InfoMO:
    """Read the Molecular orbital from the orca output."""
    _n_contracted = try_search_pattern(
        "# of contracted basis functions", file_name)
    try:
        n_contracted = _n_contracted.rsplit(maxsplit=1)[-1]  # type: ignore
    except AttributeError as ex:  # _n_contracted can be None
        raise RuntimeError(
            "Failed to extract molecular orbials from {file_name!r}") from ex

    # Parse the blocks of MOs
    start = 'MOLECULAR ORBITALS'
    end = '\n\n'
    block = parse_file(parse_section(start, end), file_name).asList()

    # split the block in lines  discarding the first line
    lines = block[0].splitlines()[1:]

    # Lines in each block
    n_contracted = int(n_contracted)
    block_lines = n_contracted + 4

    tuple_energies, tuple_coeffs = tuple(
        zip(*(read_column_orbitals(xs)
              for xs in chunked(lines, block_lines))))

    return InfoMO(np.hstack(tuple_energies), np.hstack(tuple_coeffs))


def read_column_orbitals(lines: Sequence[str]) -> tuple[NDArray[f8], NDArray[f8]]:
    """Read a set of maximum 6 columns containing the Molecular orbitals.

    the format similar to:
                        0         1         2         3         4         5
                   -19.12661  -0.94591  -0.47899  -0.35256  -0.28216   0.00756
                     2.00000   2.00000   2.00000   2.00000   2.00000   0.00000
                    --------  --------  --------  --------  --------  --------
    0H   1s        -0.002625 -0.226108 -0.331594 -0.177649  0.000002 -0.074461
    0H   2s         0.002789 -0.015914 -0.107136 -0.054883  0.000001 -0.398168
    0H   1pz       -0.001404 -0.016000 -0.008447 -0.009195  0.027232  0.002729
    0H   1px        0.002631  0.023125  0.025515 -0.010664  0.005733  0.008004
    0H   1py       -0.001684 -0.022031 -0.006086 -0.022644 -0.013760  0.008735
    """
    energies = np.array(lines[1].split(), dtype=np.float64)

    coefficients = np.array(
        [z.split()[2:] for z in lines[4:]], dtype=np.float64)

    return energies, coefficients
