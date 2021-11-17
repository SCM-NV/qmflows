"""Utilities to read cp2k out files."""

import fnmatch
import logging
import mmap
import os
import subprocess
import warnings
import re
from itertools import islice, chain
from pathlib import Path
from typing import Any, Dict, FrozenSet, Generator, Iterable, List, IO
from typing import Optional as Optional_
from typing import Sequence, Tuple, Type, TypeVar, Union, overload, Iterator

import numpy as np
from more_itertools import chunked
from pyparsing import (FollowedBy, Group, Literal, NotAny, OneOrMore, Optional,
                       SkipTo, Suppress, Word, ZeroOrMore, alphanums, alphas,
                       nums, oneOf, restOfLine, srange)
from scm.plams import Molecule, Units

from ..common import AtomBasisData, AtomBasisKey, InfoMO, MO_metadata, CP2KVersion
from ..type_hints import Literal as Literal_
from ..type_hints import PathLike, T, WarnDict, WarnMap
from ..utils import file_to_context
from ..warnings_qmflows import QMFlows_Warning
from .parser import (floatNumber, minusOrplus, natural, point,
                     try_search_pattern)
from .xyzParser import manyXYZ, tuplesXYZ_to_plams

__all__ = ['readCp2KBasis', 'read_cp2k_coefficients', 'get_cp2k_freq',
           'read_cp2k_number_of_orbitals', 'read_cp2k_xyz', 'read_cp2k_table',
           'read_cp2k_table_slc', 'get_cp2k_version', 'get_cp2k_version_run']


# Starting logger
logger = logging.getLogger(__name__)


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


def read_cp2k_coefficients(
    path_mos: "str | os.PathLike[str]",
    plams_dir: "None | str | os.PathLike[str]" = None,
) -> Union[InfoMO, Tuple[InfoMO, InfoMO]]:
    """Read the MO's from the CP2K output.

    First it reads the number of ``Orbitals`` and ``Orbital`` functions from the
    cp2k output and then read the molecular orbitals.

    Returns
    -------
        NamedTuple containing the Eigenvalues and the Coefficients
    """
    plams_dir = Path(plams_dir) if plams_dir else Path(os.getcwd())

    file_in = fnmatch.filter(os.listdir(plams_dir), '*in')[0]
    file_out = fnmatch.filter(os.listdir(plams_dir), '*out')[0]
    file_run = fnmatch.filter(os.listdir(plams_dir), '*run')[0]

    path_in = plams_dir / file_in
    path_out = plams_dir / file_out
    cp2k_version = get_cp2k_version_run(plams_dir / file_run)

    orbitals_info = read_cp2k_number_of_orbitals(path_out)
    _, range_mos = read_mos_data_input(path_in)

    try:
        # Read the range of printed MOs from the input
        printed_orbitals = range_mos[1] - range_mos[0] + 1

        return read_log_file(path_mos, printed_orbitals, orbitals_info, cp2k_version)

    except ValueError as err:
        msg = (
            f"There was a problem reading the molecular orbitals from {os.fspath(path_mos)!r},"
            "contact the developers!!"
        )
        logger.error(msg, exc_info=err)
        raise
    except TypeError as err:
        msg = (
            "There was a problem reading the ``range_mos` parameter."
            f"Its value is {range_mos}"
        )
        logger.error(msg, exc_info=err)
        raise


def read_log_file(
    path: "str | os.PathLike[str]",
    norbitals: int,
    orbitals_info: MO_metadata,
    cp2k_version: Tuple[int, int] = (0, 0),
) -> Union[InfoMO, Tuple[InfoMO, InfoMO]]:
    """
    Read the orbitals from the Log file.

    Notes
    -----
    IF IT IS A UNRESTRICTED CALCULATION, THERE ARE TWO SEPARATED SET OF MO
    FOR THE ALPHA AND BETA

    Parameters
    ----------
    path
        Path to the file containing the MO coefficients
    norbitals
        Number of MO to read
    norbital_functions
        Number of orbital functions
    cp2k_version : tuple[int, int]
        The CP2K major and minor version

    Returns
    -------
        Molecular orbitals and orbital energies
    """
    # remove Molecular orbitals coming from a restart
    move_restart_coefficients_recursively(path)

    # There is a single set of MOs
    if orbitals_info.nspinstates == 1:
        return read_coefficients(path, norbitals, orbitals_info.nOrbFuns, cp2k_version)
    else:
        path_alphas, path_betas = split_unrestricted_log_file(path)
        alphas = read_coefficients(path_alphas, norbitals, orbitals_info.nOrbFuns, cp2k_version)
        betas = read_coefficients(path_betas, norbitals - 1, orbitals_info.nOrbFuns, cp2k_version)
        return alphas, betas


def read_coefficients(
    path: PathLike,
    norbitals: int,
    norbital_functions: int,
    cp2k_version: Tuple[int, int] = (0, 0),
) -> InfoMO:
    """Read the coefficients from the plain text output.

    MO coefficients are stored in Column-major order.
    CP2K molecular orbitals output looks like:

    MO EIGENVALUES, MO OCCUPATION NUMBERS, AND SPHERICAL MO EIGENVECTORS

                              5                      6
                          -0.2590267204166110    -0.1785544120250688

                           2.0000000000000000     2.0000000000000000

     1     1  C  2s        0.0021482361354044     0.0000000235522485
     2     1  C  3s       -0.0007100065367389     0.0000000102096730
     3     1  C  3py      -0.1899052318987045     0.0000000059435027
     4     1  C  3pz       0.0000000178537720    -0.5500605729231620
     5     1  C  3px       0.3686765614894165     0.0000000228716009
     6     1  C  4py       0.0014072130025389     0.0000000019199413
     7     1  C  4pz      -0.0000000014121887     0.0293850516018881
     8     1  C  4px      -0.0028383911872079     0.0000000042372601
     9     1  C  4d-2      0.0311183981707317     0.0000000014108937
    10     1  C  4d-1      0.0000000095952723     0.0253837978837068
    11     1  C  4d0       0.0005419630026689     0.0000000391888080
    12     1  C  4d+1     -0.0000000210955114     0.0147105486663415
    13     1  C  4d+2      0.0534202997324328     0.0000000021056315
    """
    if cp2k_version >= (8, 2):
        rs: Iterable[List[str]] = _get_mos_ge_82(path)
    else:
        rs = _get_mos(path)

    # Split the list in chunks containing the orbitals info
    # in block cotaining a maximum of two columns of MOs
    chunks = chunked(rs, norbital_functions + 3)

    energies = np.empty(norbitals)
    coefficients = np.empty((norbital_functions, norbitals))

    for i, lines in enumerate(chunks):
        j = 2 * i
        es = lines[1]
        # There could be a single or double column
        css = [k[4:] for k in lines[3:]]
        # There is an odd number of MO and this is the last one
        if len(es) == 1:
            energies[-1] = float(es[0])
            coefficients[:, -1] = np.concatenate(css)
        else:
            # rearrange the coefficients
            css2 = np.transpose(css)
            energies[j: j + 2] = es
            coefficients[:, j] = css2[0]
            coefficients[:, j + 1] = css2[1]

    return InfoMO(energies, coefficients)


def _get_mos(path: PathLike) -> List[List[str]]:
    """Parse CP2k <8.2 MOs."""
    with open(path, 'r') as f:
        rs = list(filter(None, (x.rstrip('\n').split() for x in f)))
    return remove_trailing(rs[1:])


def _get_mos_ge_82(path: PathLike) -> Iterator[List[str]]:
    """Parse CP2k >=8.2 MOs."""
    with open(path, 'r') as f:
        xs = f.read().split("\n MO|")
    lineno_range = []

    # Find the begining of the MO-range
    #
    # Note that, depending on the type of CP2K executable, multiple headers may be present
    # (this seems to be a bug?)
    header = " EIGENVALUES, OCCUPATION NUMBERS, AND SPHERICAL EIGENVECTORS"
    for i, item in enumerate(xs, start=1):
        if item == header:
            lineno_range.append(i)
    if not lineno_range:
        raise ValueError("Failed to identify the start of the MO range")

    # Find the end of the MO-range
    footer_list = ["Fermi", "HOMO-LUMO", "Band"]
    for j, item in enumerate(reversed(xs)):
        if not item or any(footer in item for footer in footer_list):
            continue
        lineno_range.append(i - j)
        break
    else:
        raise ValueError("Failed to identify the end of the MO range")

    # Only read the relevant MOs
    *_, start, stop = lineno_range
    return (x.split() for x in islice(xs, start, stop) if x)


def remove_trailing(xs: List[List[str]]) -> List[List[str]]:
    """Remove the last lines of the MOs output."""
    words = {'Fermi', 'HOMO-LUMO'}
    if any(x in words for x in xs[-1]):
        xs.pop(-1)
        return remove_trailing(xs)
    else:
        return xs

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

parser_atom_label = (
    Word(srange("[A-Z]"), max=1) +
    Optional(Word(srange("[a-z]"), max=1))
)

parser_basis_name = Word(alphanums + "-") + Suppress(restOfLine)

parser_format = OneOrMore(natural + NotAny(FollowedBy(point)))

parser_key = (
    parser_atom_label.setResultsName("atom") +
    parser_basis_name.setResultsName("basisName") +
    Suppress(Literal("1"))
)

parser_basis_data = OneOrMore(floatNumber)

parser_basis = (
    parser_key +
    parser_format.setResultsName("format") +
    parser_basis_data.setResultsName("coeffs")
)

top_parser_basis = (
    OneOrMore(Suppress(comment)) + OneOrMore(
        Group(parser_basis + Suppress(Optional(OneOrMore(comment)))))
)


# ===============================<>====================================
# Parsing From File

#: A tuple with 2 elements; output of :func:`read_mos_data_input`.
Tuple2 = Tuple[Optional_[int], Optional_[Tuple[int, int]]]


def read_mos_data_input(path_input: PathLike) -> Tuple2:
    """Try to read the added_mos parameter and the range of printed MOs."""
    l1 = try_search_pattern("ADDED_MOS", path_input)
    l2 = try_search_pattern("MO_INDEX_RANGE", path_input)

    added_mos = int(l1.split()[-1]) if l1 is not None else None
    range_mos = tuple(map(int, l2.split()[1:])) if l2 is not None else None

    return added_mos, range_mos


def read_cp2k_number_of_orbitals(file_name: PathLike) -> MO_metadata:
    """Look for the line ' Number of molecular orbitals:'."""
    def fun_split(string: Optional_[str]) -> int:
        if string is not None:
            return int(string.rsplit(maxsplit=1)[-1])
        else:
            return 0

    properties = ["Number of occupied orbitals", "Number of molecular orbitals",
                  "Number of orbital functions"]

    values = [fun_split(try_search_pattern(x, file_name)) for x in properties]

    # Search for the spin states
    spin_1 = try_search_pattern("Spin 1", file_name)
    spin_2 = try_search_pattern("Spin 2", file_name)
    if all(x is not None for x in (spin_1, spin_2)):
        # There are two spin states
        values.append(2)  # It is a triplet state

    # construct the metadata inmutable object
    meta = MO_metadata(*values)

    if not all(meta):
        raise RuntimeError(f"Failed to identify orbitals in {file_name!r}")

    return meta


def move_restart_coefficients_recursively(path: "str | os.PathLike[str]") -> None:
    """Remove all the coefficients belonging to a restart."""
    while True:
        with open(path, 'r') as f:
            has_restart_coefficients = "AFTER SCF STEP -1" in f.read()

        if has_restart_coefficients:
            move_restart_coeff(path)
        else:
            break


def move_restart_coeff(path: "str | os.PathLike[str]") -> None:
    """Rename Molecular Orbital Coefficients and EigenValues."""
    root, file_name = os.path.split(path)

    # Split File into the old and new set of coefficients
    split_log_file(path, root, file_name)

    # Move the new set of coefficients to the Log file
    os.rename(Path(root, 'coeffs1'), path)

    # Remove old set of coefficients
    os.remove(Path(root, 'coeffs0'))


def split_unrestricted_log_file(path: "str | os.PathLike[str]") -> Tuple[Path, Path]:
    """Split the log file into alpha and beta molecular orbitals."""
    _root, file_name = os.path.split(path)
    split_log_file(path, _root, file_name)

    # Check that the files exists
    root = Path(_root)
    predicate = all((root / f).exists() for f in ('coeffs0', 'coeffs1'))
    if not predicate:
        msg = "There is a problem splitting the coefficients in alpha and beta components!"
        raise RuntimeError(msg)

    alpha_orbitals = Path(root, "alphas.log")
    beta_orbitals = Path(root, "betas.log")

    # Move the new set of coefficients to their corresponding names
    try:
        os.rename(Path(root, 'coeffs0'), alpha_orbitals)
        os.rename(Path(root, 'coeffs1'), beta_orbitals)
    except FileNotFoundError as ex:
        msg = "There is a problem splitting the coefficients in alpha and beta components!"
        raise RuntimeError(msg) from ex

    return alpha_orbitals, beta_orbitals


def split_log_file(path: PathLike, root: str, file_name: str) -> None:
    """Split the log file into two files with their own orbitals."""
    string = "HOMO-LUMO" if is_string_in_file("HOMO-LUMO", path) else "Fermi"
    cmd = f'csplit -f coeffs -n 1 {file_name} "/{string}/+2"'
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, cwd=root)


#: A 2-tuple; output of :func:`readCp2KBasis`.
Tuple2List = Tuple[List[AtomBasisKey], List[AtomBasisData]]


def readCp2KBasis(path: PathLike) -> Tuple2List:
    """Read the Contracted Gauss function primitives format from a text file."""
    bss = top_parser_basis.parseFile(path)
    atoms = [''.join(xs.atom[:]).lower() for xs in bss]
    names = [' '.join(xs.basisName[:]).upper() for xs in bss]
    formats = [list(map(int, xs.format[:])) for xs in bss]

    # for example 2 0 3 7 3 3 2 1 there are sum(3 3 2 1) =9 Lists
    # of Coefficients + 1 lists of exponents
    nCoeffs = [int(sum(xs[4:]) + 1) for xs in formats]
    coefficients = [list(map(float, cs.coeffs[:])) for cs in bss]
    rss = [swap_coefficients(*args) for args in zip(nCoeffs, coefficients)]
    tss = [get_head_and_tail(xs) for xs in rss]
    basisData = [AtomBasisData(xs[0], xs[1]) for xs in tss]
    basiskey = [AtomBasisKey(at, name, fmt) for at, name, fmt in zip(atoms, names, formats)]

    return (basiskey, basisData)


#: A :class:`~collections.abc.Sequence` typevar.
ST = TypeVar('ST', bound=Sequence)


@overload
def swap_coefficients(n: Literal_[1], rs: ST) -> ST:  # type: ignore
    ...


@overload
def swap_coefficients(n: int, rs: ST) -> List[ST]:
    ...


def swap_coefficients(n, rs):
    if n == 1:
        return rs
    else:
        return [rs[i::n] for i in range(n)]


def get_head_and_tail(xs: Iterable[T]) -> Tuple[T, List[T]]:
    """Return the head and tail from a list."""
    head, *tail = xs
    return (head, tail)


def get_cp2k_freq(file: Union[PathLike, IO[Any]],
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


def is_string_in_file(string: str, path: PathLike) -> bool:
    """Check if ``string`` is in file.

    .. Note::
       The search is case sensitive.
    """
    with open(path, 'r') as handler:
        s_mmap = mmap.mmap(handler.fileno(), 0, access=mmap.ACCESS_READ)
        return s_mmap.find(string.encode()) != -1


def read_cp2k_xyz(path: PathLike, dtype: Any = np.float64) -> np.ndarray:
    """Extract a 3D array from **path** with the atomic forces of all molecules.

    Requires a CP2K ``*.xyz`` file.

    """
    with open(path, 'r') as f:
        n_atom = int(next(f))
        flat_iter = chain.from_iterable(_read_cp2k_xyz(f, n_atom))
        ret = np.fromiter(flat_iter, dtype=dtype)
        ret.shape = -1, n_atom, 3  # (n_mol, n_atom, 3)
    return ret


def _read_cp2k_xyz(f: Iterable[str], n_atom: int) -> Generator[Iterator[str], None, None]:
    """Create a generator for :func:`read_cp2k_xyz`."""
    stop = 1 + n_atom
    # Account for the fact that `read_cp2k_xyz` already iterated through
    # the first element
    yield chain.from_iterable(at.split()[1:] for at in islice(f, 1, stop))
    for _ in f:
        yield chain.from_iterable(at.split()[1:] for at in islice(f, 1, stop))


def read_cp2k_table(
    path: PathLike,
    column: int,
    start: Optional_[int] = None,
    stop: Optional_[int] = None,
    step: Optional_[int] = None,
    dtype: Any = np.float64,
) -> np.ndarray:
    """Extract a 1D array from the specified **column** in **path**.

    **start**, **stop** and **step** can be used for specifiying the to-be parsed rows.

    """
    with open(path, 'r') as f:
        flat_iter = (i.split()[column] for i in islice(f, start, stop, step))
        return np.fromiter(flat_iter, dtype=dtype)


def read_cp2k_table_slc(
    path: PathLike,
    shape: Sequence[int],
    column_start: Optional_[int] = None,
    column_stop: Optional_[int] = None,
    column_step: Optional_[int] = None,
    row_start: Optional_[int] = None,
    row_stop: Optional_[int] = None,
    row_step: Optional_[int] = None,
    dtype: Any = np.float64,
) -> np.ndarray:
    """Extract an ND array of the given **shape** from **path**.

    **column_start**, **column_stop** and **column_step** can be used for
    specifiying the to-be parsed columns.
    **start**, **stop** and **step** can be used for specifiying the to-be parsed rows.

    """
    with open(path, 'r') as f:
        clm_slc = slice(column_start, column_stop, column_step)
        row_slc = islice(f, row_start, row_stop, row_step)
        flat_iter = chain.from_iterable(i.split()[clm_slc] for i in row_slc)
        return np.fromiter(flat_iter, dtype=dtype).reshape(shape)


def _get_pressure_iter(major: int, f: IO[str]) -> Generator[str, None, None]:
    """Helper function for :func:`read_cp2k_pressure`."""
    # NOTE: CP2K 8.* changed the strucure of its `.out` files,
    # hence the different prefix
    prefix1 = " MD_PAR| Pressure" if major >= 8 else " MD| Pressure"
    prefix2 = " MD| Pressure" if major >= 8 else " PRESSURE"

    # Read the initial pressure
    for i in f:
        if i.startswith(prefix1):
            yield i.split()[-1]
            break
    else:
        raise RuntimeError("Failed to identify the initial pressure")

    # Read all subsequent pressures
    for i in f:
        if i.startswith(prefix2):
            yield i.split()[-2]


def read_cp2k_pressure(
    path: PathLike,
    start: Optional_[int] = None,
    stop: Optional_[int] = None,
    step: Optional_[int] = None,
    dtype: Any = np.float64,
) -> np.ndarray:
    """Return all pressures from the passed cp2k ``.out`` file as an array."""
    major, _ = get_cp2k_version(path)

    # Read the pressures
    with open(path, 'r') as f:
        iterator = _get_pressure_iter(major, f)
        return np.fromiter(islice(iterator, start, stop, step), dtype=dtype)


def get_cp2k_version(out_file: PathLike) -> CP2KVersion:
    """Read the CP2K major and minor version from the passed .out file.

    Returns :code:`(0, 0)` if the versions cannot be identified.
    """
    with open(out_file, 'r') as f:
        for i in f:
            if i.startswith(" CP2K| version string:"):
                version_str = i.split()[-1]

                # if an error is encoutered here then we must be dealing with
                # a very old CP2K version; fall back to `major = 0` in such case
                try:
                    return CP2KVersion._make(int(i) for i in version_str.split("."))
                except ValueError:
                    pass
    warnings.warn("Failed to identify the CP2K version", QMFlows_Warning, stacklevel=2)
    return CP2KVersion(0, 0)


EXECUTABLE_PATTERN = re.compile(r"""(".+"|'.+'|\S+)\s+-i""")
VERSION_PATTERN = re.compile(r"CP2K version (\d+).(\d+)")


def get_cp2k_version_run(run_file: PathLike) -> CP2KVersion:
    """Get the CP2K version using the PLAMS .run."""
    # Extract the executable
    with open(run_file, 'r') as f:
        match = EXECUTABLE_PATTERN.search(f.read())
    if match is None:
        raise ValueError(f"Failed to extract the CP2K executable from {f.name!r}")
    executable = match.groups()[0]

    # Get the `--version` of the executable
    try:
        out = subprocess.run(
            f"{executable} --version",
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf8", shell=True
        )
    except subprocess.CalledProcessError as ex:
        raise ValueError(f"Failed to execute `{executable} --version`:\n\n{ex.stderr}") from ex

    # Parse the `--version` output
    match = VERSION_PATTERN.search(out.stdout)
    if match is None:
        raise ValueError(f"Failed to parse the `{executable} --version` output:\n\n{out.stdout}")
    return CP2KVersion._make(int(i) for i in match.groups())
