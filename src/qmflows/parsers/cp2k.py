"""Utilities to read cp2k out files."""

from __future__ import annotations

import os
import warnings
import re
from itertools import islice, chain
from collections.abc import Generator, Iterable, Iterator, Sequence
from typing import Any, IO, TYPE_CHECKING, Literal, cast

import numpy as np
from pyparsing import SkipTo, Suppress, ZeroOrMore
from scm.plams import Molecule, Units

from ..common import CP2KVersion
from ..type_hints import PathLike, WarnDict, WarnMap
from ..utils import file_to_context
from ..warnings_qmflows import QMFlows_Warning
from ._xyz import manyXYZ, tuplesXYZ_to_plams

# Re-exports
from ._cp2k_basis_parser import read_cp2k_basis
from ._cp2k_orbital_parser import read_cp2k_coefficients

if TYPE_CHECKING:
    from numpy.typing import NDArray, DTypeLike
    from numpy import float64 as f8

    #: A generator yielding 2-tuples with warning messages and warning types.
    WarnGenerator = Generator[tuple[str, type[Warning]], None, None]
    Quantity = Literal['E', 'ZPE', 'H', 'S', 'G']

__all__ = ['read_cp2k_basis', 'read_cp2k_coefficients', 'get_cp2k_freq',
           'read_cp2k_xyz', 'read_cp2k_table',
           'read_cp2k_table_slc', 'get_cp2k_version']


def read_xyz_file(file_name: PathLike) -> Molecule:
    """Read the last geometry from the output file."""
    geometries = manyXYZ(file_name)
    return tuplesXYZ_to_plams(geometries[-1])


def parse_cp2k_warnings(
    file_name: str | os.PathLike[str],
    package_warnings: WarnMap,
) -> None | WarnDict:
    """Parse All the warnings found in an output file."""
    warnings: WarnDict = {}
    for msg, named_tup in package_warnings.items():
        msg_list = named_tup.parser.parseFile(os.fspath(file_name)).asList()

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


def assign_warning(warning_type: type[Warning], msg: str,
                   msg_list: Iterable[str]) -> WarnGenerator:
    """Assign an specific Warning from the ``package_warnings`` or a generic warnings."""
    for m in msg_list:
        if msg in m:
            yield m, warning_type
        else:
            yield m, QMFlows_Warning


def get_cp2k_freq(file: PathLike | IO[Any],
                  unit: str = 'cm-1', **kwargs: Any) -> NDArray[f8]:
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
        ret = np.fromiter(iterator, dtype=np.float64, count=count)
        ret *= Units.conversion_ratio('cm-1', unit)
        return ret


QUANTITY_MAPPING = {
    ' VIB|              Electronic energy (U) [kJ/mol]:': 'E',
    ' VIB|              Zero-point correction [kJ/mol]:': 'ZPE',
    ' VIB|              Enthalpy correction (H-U) [kJ/mol]:': 'H',
    ' VIB|              Entropy [kJ/(mol K)]:': 'S',
    'VIB|              Temperature [K]:': 'T'
}

QUANTITY_SET: frozenset[Quantity] = frozenset({'E', 'ZPE', 'H', 'S', 'G'})


def get_cp2k_thermo(
    file_name: str | os.PathLike[str],
    quantity: Quantity = 'G',
    unit: str = 'kcal/mol',
) -> float:
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
    quantity = cast("Quantity", quantity.upper())
    if quantity not in QUANTITY_SET:
        raise ValueError(f"'quantity' has an invalid value ({quantity!r}); "
                         f"expected values: {tuple(QUANTITY_SET)!r}")

    parser = ZeroOrMore(Suppress(SkipTo(" VIB|              Temperature ")) + SkipTo('\n\n\n\n'))

    energy = next(iter(parser.parseFile(os.fspath(file_name))))
    energy_iter = (i.rsplit(maxsplit=1) for i in energy.splitlines() if i)
    energy_dict = {QUANTITY_MAPPING.get(k): float(v) for k, v in energy_iter}

    energy_dict['H'] += energy_dict['E']
    energy_dict['ZPE'] += energy_dict['E']
    energy_dict['G'] = energy_dict['H'] - energy_dict['T'] * energy_dict['S']

    return energy_dict[quantity] * Units.conversion_ratio('kj/mol', unit)


def read_cp2k_xyz(path: PathLike, dtype: DTypeLike = np.float64) -> NDArray[Any]:
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
    start: None | int = None,
    stop: None | int = None,
    step: None | int = None,
    dtype: DTypeLike = np.float64,
) -> NDArray[Any]:
    """Extract a 1D array from the specified **column** in **path**.

    **start**, **stop** and **step** can be used for specifiying the to-be parsed rows.

    """
    with open(path, 'r') as f:
        flat_iter = (i.split()[column] for i in islice(f, start, stop, step))
        return np.fromiter(flat_iter, dtype=dtype)


def read_cp2k_table_slc(
    path: PathLike,
    shape: Sequence[int],
    column_start: None | int = None,
    column_stop: None | int = None,
    column_step: None | int = None,
    row_start: None | int = None,
    row_stop: None | int = None,
    row_step: None | int = None,
    dtype: DTypeLike = np.float64,
) -> NDArray[Any]:
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
    start: None | int = None,
    stop: None | int = None,
    step: None | int = None,
    dtype: DTypeLike = np.float64,
) -> NDArray[Any]:
    """Return all pressures from the passed cp2k ``.out`` file as an array."""
    major, _ = get_cp2k_version(path)

    # Read the pressures
    with open(path, 'r') as f:
        iterator = _get_pressure_iter(major, f)
        return np.fromiter(islice(iterator, start, stop, step), dtype=dtype)


PATTERN = re.compile(r"CP2K version ([0-9]+)\.([0-9]+)")


def get_cp2k_version(out_file: PathLike) -> CP2KVersion:
    """Read the CP2K major and minor version from the passed .out file.

    Returns :code:`(0, 0)` if the versions cannot be identified.
    """
    with open(out_file, 'r') as f:
        for i in f:
            i = i.strip()
            if i.startswith("CP2K| version string:"):
                match = PATTERN.search(i)
                if match is None:
                    continue
                return CP2KVersion(int(match[1]), int(match[2]))

    warnings.warn(
        f"Failed extract the CP2K version from {os.fsdecode(out_file)!r}",
        QMFlows_Warning, stacklevel=2,
    )
    return CP2KVersion(0, 0)
