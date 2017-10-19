
__all__ = [
    'parse_basis_set', 'parse_hessian', 'parse_frequencies', 'parse_molecule',
    'parse_molecular_orbitals', 'parse_molecule_traj', 'parse_normal_modes']

from itertools import chain
from scm.plams import (Atom, Molecule)
from pyparsing import (alphanums, Group, OneOrMore, Word)
from .parser import (floatNumber, natural, parse_file, parse_section, skipLine,
                     skipSupress, string_array_to_molecule, try_search_pattern)
from qmflows.common import (AtomBasisData, AtomBasisKey, InfoMO)
from qmflows.utils import chunksOf
from .xyzParser import manyXYZ

import numpy as np
import pyparsing as pa

# Type hints
from typing import (List, Tuple)

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


def parse_molecular_orbitals(file_name: str) -> Tuple:
    """
    Read the Molecular orbital from the orca output
    """
    n_contracted = try_search_pattern(
        "# of contracted basis functions", file_name).split()[-1]

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
              for xs in chunksOf(lines, block_lines))))

    return InfoMO(np.hstack(tuple_energies), np.hstack(tuple_coeffs))


def read_column_orbitals(lines: List) -> Tuple:
    """
    Read a set of maximum 6 columns containing the Molecular orbitals in a
      format similar to:
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
    energies = np.array(lines[1].split(), dtype=np.float)

    coefficients = np.array(
        [z.split()[2:] for z in lines[4:]], dtype=np.float)

    return energies, coefficients


def parse_basis_set(file_name: str) -> Tuple:
    """
    Read the basis set used by Orca. It is printed by specifying the keyword:
      !printbase
    """
    parse_basis_name = skipSupress('Your calculation utilizes the basis: ') + \
        skipSupress(pa.Literal(':')) + pa.Suppress(pa.Literal(':')) + \
        pa.SkipTo('\n')
    header = skipSupress('BASIS SET IN INPUT FORMAT') + skipLine * 2
    parserElements = create_parser_element()
    parser = parse_basis_name + pa.Suppress(header) + \
        pa.Group(pa.OneOrMore(parserElements))

    basis = parse_file(parser, file_name).asList()

    return create_basis_data(basis)


def create_basis_data(basis: List) -> Tuple:
    """
    Convert the parse data into Contracted Gauss functions information
    """
    basis_name = basis[0]
    atom_keys, atom_basis = zip(*[create_CGFs_per_atom(basis_name, xs)
                                  for xs in basis[1]])

    return atom_keys, atom_basis


def create_CGFs_per_atom(basis_name: str, xs: List) -> Tuple:
    """
    Create the structure of the CGFs for each atom
    """
    atom_name = xs[0]
    formats, primitives = zip(*[((x[0], int(x[1])), x[2:]) for x in xs[1:]])

    # flatten de primitives
    primitives = list(chain(*chain(*primitives)))

    # order in coefficients and exponents the CGFs
    nPrimitives = len(primitives) // 2
    exponents, coefficients = np.array(
        primitives, dtype=np.float).reshape(nPrimitives, 2).transpose()

    atom_basis_key = AtomBasisKey(atom_name, basis_name, formats)
    atom_basis_data = AtomBasisData(exponents, coefficients)

    return atom_basis_key, atom_basis_data


def create_parser_element():
    """
    Parser to read the basis set of a given element in
    the following format:

    # Basis set for element : O
      NewGTO O
      S 5
        1    2266.1767785000     -0.0053893504
        2     340.8701019100     -0.0402347214
        3      77.3631351670     -0.1800818421
        4      21.4796449400     -0.4682885766
        5       6.6589433124     -0.4469261716
      S 1
        1       0.8097597567      1.0000000000
      S 1
        1       0.2553077223      1.0000000000
      P 3
        1      17.7215043170      0.0626302488
        2       3.8635505440      0.3333113849
        3       1.0480920883      0.7414863830
      P 1
        1       0.2764154441      1.0000000000
      D 1
        1       1.2000000000      1.0000000000
       end;
    """
    header = pa.Suppress(pa.Literal('# Basis set for element :')) + \
        pa.Word(pa.alphas) + skipLine

    parseCGF = pa.Group(
        pa.Word(pa.alphas, exact=1)  + natural +
        pa.Group(pa.OneOrMore(pa.Suppress(natural) + floatNumber * 2)))

    return pa.Group(header + pa.OneOrMore(parseCGF) + pa.Suppress(pa.Literal('end;')))
