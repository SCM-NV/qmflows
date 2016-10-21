
from plams import (Atom, Molecule)
from pyparsing import (alphanums, Group, OneOrMore, Word)
from .parser import (floatNumber, parse_file, parse_section, skipLine,
                     skipSupress, string_array_to_molecule)
from qmworks.utils import chunksOf
from .xyzParser import manyXYZ

import numpy as np

# Type hints
Vector = np.ndarray
Matrix = np.ndarray

vectorize_float = np.vectorize(float)


def parse_molecule(file_name, mol=None):
    """
    Parse The Cartesian coordinates from the output file.
    """
    header = "CARTESIAN COORDINATES (ANGSTROEM)"
    p1 = skipSupress(header) + skipLine * 2
    parse_atoms = Group(Word(alphanums) + floatNumber * 3)
    parse_mol = p1 + Group(OneOrMore(parse_atoms))

    parse_many_mol = OneOrMore(parse_mol)

    return string_array_to_molecule(parse_many_mol, file_name, mol=mol)


def parse_molecule_traj(file_traj):
    """
    Read Molecules from the job_name.traj file.
    """
    mols = manyXYZ(file_traj)
    # Last geometry corresponds to the optimized structure
    opt_mol = mols[-1]

    plams_mol = Molecule()
    for at in opt_mol:
        symb = at.symbol
        cs = at.xyz
        plams_mol.add_atom(Atom(symbol=symb, coords=tuple(cs)))

    return plams_mol


def parse_hessian(file_hess: str) -> Matrix:
    """
    Read the hessian matrix in cartesian coordinates from the job_name.hess file.
    :returns: Numpy array
    """
    start = '$hessian'
    return read_blocks_from_file(start, '\n\n', file_hess)


def parse_normal_modes(file_hess: str) -> Matrix:
    """
    Returns the normal modes from the job_name.hess file
    """
    start = '$normal_modes'
    return read_blocks_from_file(start, '\n\n', file_hess)


def parse_frequencies(file_hess: str) -> Matrix:
    """
    Parse the vibrational frequencies from the job_name.hess file.
    """
    p = parse_section('$vibrational_frequencies', '\n\n')
    lines = parse_file(p, file_hess)[0].splitlines()

    return vectorize_float([x.split()[-1] for x in lines[1:]])


def read_blocks_from_file(start: str, end: str, file_name: str) -> Matrix:
    """
    Read a matrix printed in block format

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
              for block in chunksOf(lines, number_of_basis + 1)]

    return np.concatenate(blocks, axis=1)


def read_block(lines):
    """
    Read a block containing the values of the block matrix in
    a format similar to:

              0          1          2          3          4          5
      0   0.752994  0.078644   0.000000  -0.134007   0.156950  -0.000000
      1   0.078857  0.848881  -0.000000   0.162989  -0.312155   0.000000
      2   0.000000 -0.000000   0.074907  -0.000000   0.000000  -0.026031
      3  -0.133981  0.162928  -0.000000   0.133434  -0.164602   0.000000

    Read the matrix skiping the header and the first integer index

    :returns: Numpy array
    """
    return np.stack(map(lambda x: vectorize_float(x.split()[1:]), lines[1:]))
