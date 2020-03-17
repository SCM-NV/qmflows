"""Utilities to read cp2k out files."""

__all__ = ['readCp2KBasis', 'read_cp2k_coefficients',
           'read_cp2k_number_of_orbitals']

import fnmatch
import os
import subprocess
from collections import namedtuple
from itertools import islice
from typing import Type, Generator, Iterable, Tuple
from typing import Optional as Optional_

import numpy as np
from more_itertools import chunked
from pyparsing import (
    FollowedBy, Group, Literal, NotAny, OneOrMore, Optional,
    Suppress, Word, alphanums, alphas, nums, oneOf, restOfLine, srange
)

from .parser import floatNumber, minusOrplus, natural, point, try_search_pattern
from .xyzParser import manyXYZ, tuplesXYZ_to_plams
from ..common import AtomBasisData, AtomBasisKey, InfoMO
from ..warnings_qmflows import QMFlows_Warning
from ..type_hints import WarnMap, WarnDict, PathLike

# =========================<>=============================
MO_metadata = namedtuple("MO_metadada", ("nOccupied", "nOrbitals", "nOrbFuns"))

# Molecular Orbitals Parsing
# MO EIGENVALUES, MO OCCUPATION NUMBERS, AND SPHERICAL MO EIGENVECTORS

# 1                      2
#                          -0.9857682741370732    -0.9831467097855797

#                           2.0000000000000000     2.0000000000000000

#     1     1 cd  2s       -0.0015026981889089    -0.0103313715516893
#     2     1 cd  3s       -0.0005376142747880    -0.0041729598190025
#     3     1 cd  3py      -0.0013790317507575     0.0132729535025288
#     4     1 cd  3pz      -0.0015557487597731    -0.0005486094359245
#     5     1 cd  3px      -0.0013339995106232    -0.0100914249163043
#     6     1 cd  4py      -0.0003884918433452     0.0046068283721132


def read_xyz_file(file_name: PathLike):
    """Read the last geometry from the output file."""
    geometries = manyXYZ(file_name)
    return tuplesXYZ_to_plams(geometries[-1])


def parse_cp2k_warnings(file_name: PathLike,
                        package_warnings: WarnMap) -> Optional_[WarnDict]:
    """Parse All the warnings found in an output file."""
    warnings = {}
    for msg, named_tup in package_warnings.items():
        msg_list = named_tup.parser.parseFile(file_name).asList()

        # Search for warnings that match the ones provided by the user
        iterator = assign_warning(named_tup.warn_type, msg, msg_list)

        # Apply post processing to the exception message
        for msg_ret, warn_type in iterator:
            key = named_tup.func(msg_ret)
            if key is not None:
                v = warnings.get(key, QMFlows_Warning)
                if v is QMFlows_Warning:
                    warnings[key] = warn_type

    return warnings or None


def assign_warning(warning_type: Type[Warning], msg: str,
                   msg_list: Iterable[str]
                   ) -> Generator[Tuple[str, Type[Warning]], None, None]:
    """Assign an specific Warning from the ``package_warnings`` or a generic warnings."""
    for m in msg_list:
        if msg in m:
            yield m, warning_type
        else:
            yield m, QMFlows_Warning


def read_cp2k_coefficients(path_mos, plams_dir=None):
    """Read the MO's from the CP2K output.

    First it reads the number of ``Orbitals`` and ``Orbital`` functions from the
    cp2k output and then read the molecular orbitals.

    :returns: NamedTuple containing the Eigenvalues and the Coefficients
    """
    file_out = fnmatch.filter(os.listdir(plams_dir), '*out')[0]
    file_in = fnmatch.filter(os.listdir(plams_dir), '*in')[0]
    path_in, path_out = [os.path.join(plams_dir, x)
                         for x in [file_in, file_out]]
    orbitals_info = read_cp2k_number_of_orbitals(path_out)
    added_mos, range_mos = read_mos_data_input(path_in)

    # Read the range of printed MOs from the input
    if range_mos is not None:
        printed_orbitals = range_mos[1] - range_mos[0] + 1
    # Otherwise read the added_mos parameter
    elif added_mos is not None:
        printed_orbitals = orbitals_info.added_mos * 2
    # Otherwise read the occupied orbitals
    else:
        printed_orbitals = orbitals_info.nOccupied

    return readCp2KCoeff(path_mos, printed_orbitals, orbitals_info.nOrbFuns)


def readCp2KCoeff(path, nOrbitals: int, nOrbFuns: int):
    """Read the coefficients from the plain text output.

    MO coefficients are stored in Column-major order.

    :parameter path: Path to the file containing the MO coefficients
    :type path: String
    :parameter nOrbitals: Number of MO to read
    :param nOrbFuns: Number of orbital functions
    :returns: Molecular orbitals and orbital energies
    """
    def remove_trailing(xs):
        "Remove the last lines of the MOs output."
        words = ['Fermi', 'HOMO-LUMO']
        if any([x in words for x in xs[-1]]):
            xs.pop(-1)
            return remove_trailing(xs)
        else:
            return xs

    # Check if the Molecular orbitals came from a restart
    with open(path, 'r') as f:
        xs = list(islice(f, 4))
    if "AFTER SCF STEP -1" in ''.join(xs):
        move_restart_coeff(path)

    # Open the MO file
    with open(path, 'r') as f:
        xss = f.readlines()

    # remove empty lines and comments
    rs = list(filter(None, map(lambda x: x.split(), xss)))
    rs = remove_trailing(rs[1:])  # remove header and trail comments

    # Split the list in chunks containing the orbitals info
    # in block cotaining a maximum of two columns of MOs
    chunks = chunked(rs, nOrbFuns + 3)

    eigenVals = np.empty(nOrbitals)
    coefficients = np.empty((nOrbFuns, nOrbitals))

    convert_to_float = np.vectorize(float)
    for i, xs in enumerate(chunks):
        j = 2 * i
        es = xs[1]
        css = [l[4:] for l in xs[3:]]
        # There is an odd number of MO and this is the last one
        if len(es) == 1:
            eigenVals[-1] = float(es[0])
            coefficients[:, -1] = convert_to_float(np.concatenate(css))
        else:
            # rearrange the coeff
            css = np.transpose(convert_to_float(css))
            eigenVals[j: j + 2] = es
            coefficients[:, j] = css[0]
            coefficients[:, j + 1] = css[1]

    return InfoMO(eigenVals, coefficients)

# =====================> Orbital Parsers <===================


xyz = oneOf(['x', 'y', 'z'])

orbS = Literal("s")

orbP = Literal("p") + xyz

orbD = Literal("d") + (Literal('0') |
                       (minusOrplus + Word(srange("[1-2]"), max=1)) |
                       (xyz + oneOf(['2', '3', 'y', 'z'])))

orbF = Literal("f") + (Literal('0') |
                       (minusOrplus + Word(srange("[1-3]"), max=1)) |
                       (xyz + oneOf(['2', '3', 'y', 'z']) +
                        Optional(oneOf(['2', 'y', 'z']))))

orbitals = Word(nums, max=1) + (orbS | orbP | orbD | orbF)

# Orbital Information:"        12     1 cd  4d+1"
orbInfo = natural * 2 + Word(alphas, max=2) + orbitals


def funCoefficients(x):
    """Parser Coeffcients."""
    fun = OneOrMore(Suppress(orbInfo) + floatNumber * x)
    return fun.setResultsName("coeffs")


def funOrbNumber(x):
    """Parse Orbital Occupation Number. There is min 1 max 4."""
    return natural * x


# ====================> Basis File <==========================
comment = Literal("#") + restOfLine

parseAtomLabel = (Word(srange("[A-Z]"), max=1) +
                  Optional(Word(srange("[a-z]"), max=1)))

parserBasisName = Word(alphanums + "-") + Suppress(restOfLine)

parserFormat = OneOrMore(natural + NotAny(FollowedBy(point)))

parserKey = parseAtomLabel.setResultsName("atom") + \
    parserBasisName.setResultsName("basisName") + \
    Suppress(Literal("1"))

parserBasisData = OneOrMore(floatNumber)

parserBasis = parserKey + parserFormat.setResultsName("format") + \
    parserBasisData.setResultsName("coeffs")


topParseBasis = OneOrMore(Suppress(comment)) + \
    OneOrMore(Group(parserBasis + Suppress(Optional(OneOrMore(comment)))))


# ===============================<>====================================
# Parsing From File

def read_mos_data_input(path_input):
    """Try to read the added_mos parameter and the range of printed MOs."""
    properties = ["ADDED_MOS", "MO_INDEX_RANGE"]
    l1, l2 = [try_search_pattern(x, path_input) for x in properties]
    added_mos = l1.split()[-1] if l1 is not None else None
    range_mos = list(map(int, l2.split()[1:])) if l1 is not None else None

    return added_mos, range_mos


def read_cp2k_number_of_orbitals(file_name):
    """Look for the line ' Number of molecular orbitals:'."""
    def fun_split(l):
        return l.split()[-1]

    properties = ["Number of occupied orbitals", "Number of molecular orbitals",
                  "Number of orbital functions"]

    xs = [fun_split(try_search_pattern(x, file_name)) for x in properties]

    return MO_metadata(*[int(x) for x in xs])


def move_restart_coeff(path):
    """Rename Molecular Orbital Coefficients and EigenValues."""
    root, file_name = os.path.split(path)
    # Current work directory
    cwd = os.path.realpath('.')
    # Change directory
    os.chdir(root)
    # Split File into the old and new set of coefficients
    cmd = 'csplit -f coeffs -n 1 {} "/HOMO-LUMO/+2"'.format(file_name)
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)
    # Move the new set of coefficients to the Log file
    os.rename('coeffs1', file_name)
    # Remove old set of coefficients
    os.remove('coeffs0')
    # Return to CWD
    os.chdir(cwd)


def readCp2KBasis(path: str) -> tuple:
    """Read the Contracted Gauss function primitives format from a text file."""
    bss = topParseBasis.parseFile(path)
    atoms = [''.join(xs.atom[:]).lower() for xs in bss]
    names = [' '.join(xs.basisName[:]).upper() for xs in bss]
    formats = [list(map(int, xs.format[:])) for xs in bss]
    # for example 2 0 3 7 3 3 2 1 there are sum(3 3 2 1) =9 Lists
    # of Coefficients + 1 lists of exponents
    nCoeffs = [int(sum(xs[4:]) + 1) for xs in formats]
    coefficients = [list(map(float, cs.coeffs[:])) for cs in bss]
    rss = [swapCoeff(*args) for args in zip(nCoeffs, coefficients)]
    tss = [headTail(xs) for xs in rss]
    basisData = [AtomBasisData(xs[0], xs[1]) for xs in tss]
    basiskey = [AtomBasisKey(*rs) for rs in zip(atoms, names, formats)]

    return (basiskey, basisData)


def swapCoeff(n, rs):
    if n == 1:
        return rs
    else:
        return [rs[i::n] for i in range(n)]


def headTail(xs):
    """Return the head and tail from a list."""
    it = iter(xs)
    head = next(it)
    tail = list(it)
    return (head, tail)
