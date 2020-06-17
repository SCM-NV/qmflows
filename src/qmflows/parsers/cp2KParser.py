"""Utilities to read cp2k out files."""

__all__ = ['readCp2KBasis', 'read_cp2k_coefficients', 'get_cp2k_freq',
           'read_cp2k_number_of_orbitals']

import os
import fnmatch
import subprocess
from io import TextIOBase
from pathlib import Path
from itertools import islice
from typing import (Union, Any, Type, Generator, Iterable, Tuple,
                    List, Sequence, TypeVar, overload, FrozenSet, Dict)
from typing import Optional as Optional_

import numpy as np
from scm.plams import Units, Molecule
from more_itertools import chunked
from pyparsing import (
    FollowedBy, Group, Literal, NotAny, OneOrMore, Optional,
    Suppress, Word, alphanums, alphas, nums, oneOf, restOfLine, srange,
    ZeroOrMore, SkipTo
)

from .parser import floatNumber, minusOrplus, natural, point, try_search_pattern
from .xyzParser import manyXYZ, tuplesXYZ_to_plams
from ..utils import file_to_context
from ..common import AtomBasisData, AtomBasisKey, InfoMO, MO_metadata
from ..warnings_qmflows import QMFlows_Warning
from ..type_hints import WarnMap, WarnDict, PathLike, T
from ..type_hints import Literal as Literal_


def read_xyz_file(file_name: PathLike) -> Molecule:
    """Read the last geometry from the output file."""
    geometries = manyXYZ(file_name)
    return tuplesXYZ_to_plams(geometries[-1])


def parse_cp2k_warnings(file_name: PathLike,
                        package_warnings: WarnMap) -> Optional_[WarnDict]:
    """Parse All the warnings found in an output file."""
    warnings: WarnDict = {}
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


#: A generator yielding 2-tuples with warning messages and warning types.
WarnGenerator = Generator[Tuple[str, Type[Warning]], None, None]


def assign_warning(warning_type: Type[Warning], msg: str,
                   msg_list: Iterable[str]) -> WarnGenerator:
    """Assign an specific Warning from the ``package_warnings`` or a generic warnings."""
    for m in msg_list:
        if msg in m:
            yield m, warning_type
        else:
            yield m, QMFlows_Warning


def read_cp2k_coefficients(path_mos: PathLike,
                           plams_dir: Optional_[PathLike] = None) -> InfoMO:
    """Read the MO's from the CP2K output.

    First it reads the number of ``Orbitals`` and ``Orbital`` functions from the
    cp2k output and then read the molecular orbitals.

    :returns: NamedTuple containing the Eigenvalues and the Coefficients
    """
    plams_dir = Path(plams_dir) if plams_dir else Path(os.getcwd())

    file_in = fnmatch.filter(os.listdir(plams_dir), '*in')[0]
    file_out = fnmatch.filter(os.listdir(plams_dir), '*out')[0]

    path_in = plams_dir / file_in
    path_out = plams_dir / file_out

    orbitals_info = read_cp2k_number_of_orbitals(path_out)
    added_mos, range_mos = read_mos_data_input(path_in)

    # Read the range of printed MOs from the input
    printed_orbitals = range_mos[1] - range_mos[0] + 1

    return readCp2KCoeff(path_mos, printed_orbitals, orbitals_info.nOrbFuns)


def readCp2KCoeff(path: PathLike, nOrbitals: int, nOrbFuns: int) -> InfoMO:
    """Read the coefficients from the plain text output.

    MO coefficients are stored in Column-major order.

    :parameter path: Path to the file containing the MO coefficients
    :type path: String
    :parameter nOrbitals: Number of MO to read
    :param nOrbFuns: Number of orbital functions
    :returns: Molecular orbitals and orbital energies
    """
    def remove_trailing(xs):
        """Remove the last lines of the MOs output."""
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
        css = [k[4:] for k in xs[3:]]
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

# ====================> Basis File <==========================
comment = Literal("#") + restOfLine

parseAtomLabel = (
    Word(srange("[A-Z]"), max=1) +
    Optional(Word(srange("[a-z]"), max=1))
)

parserBasisName = Word(alphanums + "-") + Suppress(restOfLine)

parserFormat = OneOrMore(natural + NotAny(FollowedBy(point)))

parserKey = (
    parseAtomLabel.setResultsName("atom") +
    parserBasisName.setResultsName("basisName") +
    Suppress(Literal("1"))
)

parserBasisData = OneOrMore(floatNumber)

parserBasis = (
    parserKey +
    parserFormat.setResultsName("format") +
    parserBasisData.setResultsName("coeffs")
)

topParseBasis = (
    OneOrMore(Suppress(comment)) +
    OneOrMore(Group(parserBasis +
    Suppress(Optional(OneOrMore(comment)))))
)


# ===============================<>====================================
# Parsing From File

#: A tuple with 2 elements; output of :func:`read_mos_data_input`.
Tuple2 = Tuple[Optional_[str], Optional_[List[int]]]


def read_mos_data_input(path_input: PathLike) -> Tuple2:
    """Try to read the added_mos parameter and the range of printed MOs."""
    l1 = try_search_pattern("ADDED_MOS", path_input)
    l2 = try_search_pattern("MO_INDEX_RANGE", path_input)

    added_mos = l1.split()[-1] if l1 is not None else None
    range_mos = list(map(int, l2.split()[1:])) if l2 is not None else None

    return added_mos, range_mos


def read_cp2k_number_of_orbitals(file_name: PathLike) -> MO_metadata:
    """Look for the line ' Number of molecular orbitals:'."""
    def fun_split(string: str) -> str:
        return string.split()[-1]

    properties = ["Number of occupied orbitals", "Number of molecular orbitals",
                  "Number of orbital functions"]

    try:
        xs = [fun_split(try_search_pattern(x, file_name)) for x in properties]  # type: ignore
    except AttributeError as ex:  # try_search_pattern() can return None
        raise RuntimeError(f"Failed to identify orbitals in {file_name!r}") from ex

    return MO_metadata(*[int(x) for x in xs])


def move_restart_coeff(path: PathLike) -> None:
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


#: A 2-tuple; output of :func:`readCp2KBasis`.
Tuple2List = Tuple[List[AtomBasisKey], List[AtomBasisData]]


def readCp2KBasis(path: PathLike) -> Tuple2List:
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
    basiskey = [AtomBasisKey(at, name, fmt) for at, name, fmt in zip(atoms, names, formats)]

    return (basiskey, basisData)


#: A :class:`~collections.abc.Sequence` typevar.
ST = TypeVar('ST', bound=Sequence)


@overload
def swapCoeff(n: Literal_[1], rs: ST) -> ST:  # type: ignore
    ...
@overload
def swapCoeff(n: int, rs: ST) -> List[ST]:
    ...
def swapCoeff(n, rs):
    if n == 1:
        return rs
    else:
        return [rs[i::n] for i in range(n)]


def headTail(xs: Iterable[T]) -> Tuple[T, List[T]]:
    """Return the head and tail from a list."""
    head, *tail = xs
    return (head, tail)


def get_cp2k_freq(file: Union[PathLike, TextIOBase],
                  unit: str = 'cm-1', **kwargs: Any) -> np.ndarray:
    r"""Extract vibrational frequencies from *file*, a CP2K .mol file in the Molden format.

    Paramters
    ---------
    file : :class:`str`, :class:`bytes`, :class:`os.PathLike` or :class:`io.IOBase`
        A `path- <https://docs.python.org/3/glossary.html#term-path-like-object>`_ or
        `file-like <https://docs.python.org/3/glossary.html#term-file-object>`_ object
        pointing to the CP2K .mol file.
        Note that passed file-like objects should return strings (not bytes) upon iteration;
        consider wrapping *file* in :func:`codecs.iterdecode` if its iteration will yield bytes.

    unit : :class:`str`
        The output unit of the vibrational frequencies.
        See :class:`plams.Units<scm.plams.tools.units.Units>` for more details.

    /**kwargs : :data:`~typing.Any`
        Further keyword arguments for :func:`open`.
        Only relevant if *file* is a path-like object.

    Returns
    -------
    :class:`numpy.ndarray` [:class:`float`], shape :math:`(n,)`
        A 1D array of length :math:`n` containing the vibrational frequencies
        extracted from *file*.

    """
    context_manager = file_to_context(file, **kwargs)

    with context_manager as f:
        item = next(f)
        if not isinstance(item, str):
            raise TypeError(f"Iteration through {f!r} should yield strings; "
                            f"observed type: {item.__class__.__name__!r}")

        # Find the start of the [Atoms] block
        elif '[Atoms]' not in item:
            for item in f:
                if '[Atoms]' in item:
                    break
            else:
                raise ValueError(f"failed to identify the '[Atoms]' substring in {f!r}")

        # Find the end of the [Atoms] block, i.e. the start of the [FREQ] block
        for atom_count, item in enumerate(f):
            if '[FREQ]' in item:
                break
        else:
            raise ValueError(f"failed to identify the '[FREQ]' substring in {f!r}")

        # Identify the vibrational degrees of freedom
        if atom_count == 0:
            raise ValueError(f"failed to identify any atoms in the '[Atoms]' section of {f!r}")
        elif atom_count <= 2:
            count = atom_count - 1
        else:
            count = 3 * atom_count - 6

        # Gather and return the frequencies
        iterator = islice(f, 0, count)
        ret = np.fromiter(iterator, dtype=float, count=count)
        ret *= Units.conversion_ratio('cm-1', unit)
        return ret


QUANTITY_MAPPING: Dict[str, str] = {
    ' VIB|              Electronic energy (U) [kJ/mol]:': 'E',
    ' VIB|              Zero-point correction [kJ/mol]:': 'ZPE',
    ' VIB|              Enthalpy correction (H-U) [kJ/mol]:': 'H',
    ' VIB|              Entropy [kJ/(mol K)]:': 'S',
    'VIB|              Temperature [K]:': 'T'
}

QUANTITY_SET: FrozenSet[str] = frozenset({'E', 'ZPE', 'H', 'S', 'G'})

Quantity = Literal_['E', 'ZPE', 'H', 'S', 'G']


def get_cp2k_thermo(file_name: PathLike, quantity: Quantity = 'G',
                    unit: str = 'kcal/mol') -> float:
    """Return thermochemical properties as extracted from a CP2K .out file.

    Note
    ----
    Note that CP2K can under certain circumstances report entropies and, by extension,
    Gibbs free energies as :data:`~math.nan`.

    Parameters
    ----------
    file_name : :class:`str`, :class:`bytes` or :class:`os.PathLike`
        A path-like object pointing to the CP2K .out file.

    quantity : :class:`str`
        The to-be returned quantity.
        Accepted values are ``"E"``, ``"ZPE"``, ``"H"``, ``"S"`` and ``"G"``.

    unit : :class:`str`
        The unit of the to-be returned *quantity*.
        See :class:`plams.Units<scm.plams.tools.units.Units>` for more details.

    Returns
    -------
    :class:`float`
        A user-specified *quantity* expressed in *unit* as extracted from *file_name*.

    """
    quantity = quantity.upper()
    if quantity not in QUANTITY_SET:
        raise ValueError(f"'quantity' has an invalid value ({quantity!r}); "
                         f"expected values: {tuple(QUANTITY_SET)!r}")

    parser = ZeroOrMore(Suppress(SkipTo(" VIB|              Temperature ")) + SkipTo('\n\n\n\n'))

    energy = next(iter(parser.parseFile(file_name)))
    energy_iter = (i.rsplit(maxsplit=1) for i in energy.splitlines() if i)
    energy_dict = {QUANTITY_MAPPING.get(k): float(v) for k, v in energy_iter}

    energy_dict['H'] += energy_dict['E']
    energy_dict['ZPE'] += energy_dict['E']
    energy_dict['G'] = energy_dict['H'] - energy_dict['T'] * energy_dict['S']

    return energy_dict[quantity] * Units.conversion_ratio('kj/mol', unit)
