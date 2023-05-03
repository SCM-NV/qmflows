"""A module with utility function for CP2K.

Index
-----
.. currentmodule:: qmflows.cp2k_utils
.. autosummary::
    set_prm
    map_psf_atoms
    prm_to_df
    CP2K_KEYS_ALIAS

API
---
.. autofunction:: set_prm
.. autofunction:: map_psf_atoms
.. autofunction:: prm_to_df
.. autodata:: CP2K_KEYS_ALIAS
    :annotation: : dict[str, tuple[str, ...]]

    .. code:: python

{cp2k_keys_alias}

"""

from __future__ import annotations

import copy
import textwrap
from functools import singledispatch
from itertools import repeat, islice
from collections.abc import MutableMapping, Sequence, Iterable, Iterator
from typing import NoReturn, Any, overload, cast, IO, TYPE_CHECKING

import numpy as np
import pandas as pd
from pandas.core.generic import NDFrame

from scm import plams

from ._settings import Settings
from .utils import to_runtime_error, file_to_context
from .type_hints import MappingScalar, MappingSequence, PathLike

if TYPE_CHECKING:
    PrmMapping = (
        MappingScalar
        | Sequence[MappingScalar]
        | MappingSequence
        | Sequence[MappingSequence]
    )

__all__ = ['set_prm', 'prm_to_df', 'map_psf_atoms', 'CP2K_KEYS_ALIAS']

_BASE_PATH = ('specific', 'cp2k', 'force_eval', 'mm', 'forcefield')

#: A dictionary mapping ``key_path`` aliases to the actual keys.
CP2K_KEYS_ALIAS: dict[str, tuple[str, ...]] = {
    'bond': _BASE_PATH + ('bond',),
    'bend': _BASE_PATH + ('bend',),
    'opbend': _BASE_PATH + ('opbend',),
    'ub': _BASE_PATH + ('bend', 'ub'),
    'torsion': _BASE_PATH + ('torsion',),
    'improper': _BASE_PATH + ('improper',),
    'charge': _BASE_PATH + ('charge',),
    'dipole': _BASE_PATH + ('dipole',),
    'quadrupole': _BASE_PATH + ('quadrupole',),

    'lennard-jones': _BASE_PATH + ('nonbonded', 'lennard-jones'),
    'lennard_jones': _BASE_PATH + ('nonbonded', 'lennard-jones'),
    'allegro': _BASE_PATH + ('nonbonded', 'allegro'),
    'bmhft': _BASE_PATH + ('nonbonded', 'bmhft'),
    'bmhftd': _BASE_PATH + ('nonbonded', 'bmhftd'),
    'buck4ranges': _BASE_PATH + ('nonbonded', 'buck4ranges'),
    'buckmorse': _BASE_PATH + ('nonbonded', 'buckmorse'),
    'eam': _BASE_PATH + ('nonbonded', 'eam'),
    'gal19': _BASE_PATH + ('nonbonded', 'gal19'),
    'gal21': _BASE_PATH + ('nonbonded', 'gal21'),
    'genpot': _BASE_PATH + ('nonbonded', 'genpot'),
    'goodwin': _BASE_PATH + ('nonbonded', 'goodwin'),
    'ipbv': _BASE_PATH + ('nonbonded', 'ipbv'),
    'nequip': _BASE_PATH + ('nonbonded', 'nequip'),
    'quip': _BASE_PATH + ('nonbonded', 'quip'),
    'siepmann': _BASE_PATH + ('nonbonded', 'siepmann'),
    'tabpot': _BASE_PATH + ('nonbonded', 'tabpot'),
    'tersoff': _BASE_PATH + ('nonbonded', 'tersoff'),
    'williams': _BASE_PATH + ('nonbonded', 'williams'),

    'lennard-jones14': _BASE_PATH + ('nonbonded14', 'lennard-jones'),
    'lennard_jones14': _BASE_PATH + ('nonbonded14', 'lennard-jones'),
    'genpot14': _BASE_PATH + ('nonbonded14', 'genpot'),
    'goodwin14': _BASE_PATH + ('nonbonded14', 'goodwin'),
    'williams14': _BASE_PATH + ('nonbonded14', 'williams'),
}
del _BASE_PATH


class LengthError(ValueError):
    """A :exc:`ValueError` subclass for exceptions caused by incorrect lengths of a :class:`~collections.abc.Mapping` or :class:`~collections.abc.Sequence`."""  # noqa: E501


@to_runtime_error
def _map_psf_atoms(
    settings: None,
    key: str,
    value: PathLike | IO[Any],
    mol: None, **kwargs: Any,
) -> dict[str, str]:
    """A small wrapper around :func:`map_psf_atoms`."""
    return map_psf_atoms(value, **kwargs)


def map_psf_atoms(file: PathLike | IO[Any], **kwargs: Any) -> dict[str, str]:
    r"""Take a .psf file and construct a :class:`dict` mapping atom types to atom names.

    Examples
    --------
    .. code:: python

        >>> from io import StringIO
        >>> from qmflows.cp2k_utils import map_psf_atoms

        >>> file = StringIO('''
        ... PSF EXT
        ...
        ...         10 !NATOM
        ...          1 MOL1     1        LIG      C        C331   -0.272182       12.010600        0
        ...          2 MOL1     1        LIG      C        C321   -0.282182       12.010600        0
        ...          3 MOL1     1        LIG      C        C2O3    0.134065       12.010600        0
        ...          4 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
        ...          5 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
        ...          6 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
        ...          7 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
        ...          8 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
        ...          9 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
        ...         10 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
        ... ''')

        >>> atom_map = map_psf_atoms(file)
        >>> print(atom_map)
        {'C331': 'C', 'C321': 'C', 'C2O3': 'C', 'O2D2': 'O', 'HGA2': 'H', 'HGA3': 'H'}


    Parameters
    ----------
    file : :class:`str`, :class:`bytes`, :class:`os.PathLike` or :class:`io.TextIOBase`
        A `path- <https://docs.python.org/3/glossary.html#term-path-like-object>`_ or
        `file-like <https://docs.python.org/3/glossary.html#term-file-object>`_ object
        pointing to the .psf file.
        Note that passed file-like objects should return strings (not bytes) upon iteration;
        consider wrapping *file* in :func:`codecs.iterdecode` if its iteration will yield bytes.

    /**kwargs : :data:`~typing.Any`
        Further keyword arguments for :func:`open`.
        Only relevant when *file* is a path-like object.

    Returns
    -------
    :class:`dict` [:class:`str`, :class:`str`]
        A dictionary mapping atom types to atom names.
        Atom types/names are extracted from the passed .psf file.

    """  # noqa: E501
    context_manager = file_to_context(file, **kwargs)

    with context_manager as f:
        # A quick check to see *f* is not opened in bytes mode or something similar
        i = next(f)
        if not isinstance(i, str):
            raise TypeError(f"Iteration through {f!r} should yield strings; "
                            f"observed type: {i.__class__.__name__!r}")
        elif '!NATOM' not in i:
            for i in f:
                if '!NATOM' in i:  # Find the !NATOM block
                    break
            else:
                raise ValueError(f"failed to identify the '!NATOM' substring in {f!r}")

        # Identify rows 4 & 5 which contain, respectivelly,
        # the `atom name` and `atom type` blocks
        atom_count = int(i.split()[0])
        iterator = (j.split()[4:6] for j in islice(f, 0, atom_count))
        try:
            return {k: v for v, k in iterator}
        except ValueError as ex:
            raise ValueError("Failed to identify the 'atom name'- and "
                             f"'atom type'-containing rows in {f!r};\n{ex}") from ex


@to_runtime_error
def set_prm(settings: Settings, key: str | tuple[str, ...],
            value: PrmMapping, mol: None | plams.Molecule) -> None:
    """Assign a set of forcefield parameters to *settings* as specific keys.

    Examples
    --------
    *value* should be a dictionary whose values:

    * Entirelly consist of scalars.
    * Entirelly consist of Sequences.

    .. code:: python

        >>> lennard_jones = {  # Example 1
        ...     'param': ('epsilon', 'sigma'),
        ...     'unit': ('kcalmol', 'angstrom'),  # An optional key
        ...     'Cs': (1, 1),
        ...     'Cd': (2, 2),
        ...     'O': (3, 3),
        ...     'H': (4, 4)
        ... }

        >>> lennard_jones = [  # Example 3
        ... {'param': 'epsilon',
        ...  'unit': 'kcalmol',  # An optional key
        ...  'Cs': 1,
        ...  'Cd': 2,
        ...  'O': 3,
        ...  'H': 4},
        ... {'param': 'sigma',
        ...  'unit': 'angstrom',  # An optional key
        ...  'Cs': 1,
        ...  'Cd': 2,
        ...  'O': 3,
        ...  'H': 4}
        ... ]


    Warning
    -------
    Scalars and sequences **cannot** be freely mixed in the values of *value*;
    it should be one or the other.

    Parameters
    ----------
    settings : :class:`qmflows.Settings`
        The input CP2K settings.

    key : :class:`str` or :class:`tuple` [:class:`str`, ...]
        A path of CP2K keys or an alias for a pre-defined path (see :data:`CP2K_KEYS_ALIAS`).

    value : :class:`~collections.abc.MutableMapping` [:class:`str`, ``T`` or :class:`~collections.abc.Sequence` [``T``]]
        A dictionary with the to-be added parameters.
        Scalars and sequences **cannot** be freely mixed in the dictionary values;
        it should be one or the other.

        See the Examples section for more details.

    mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`, optional
        A dummy argument in order to ensure signature compatiblity.

    Raises
    ------
    :exc:`RuntimeError`
        Raised when issues are encountered parsing either *key* or *value*.

    See Also
    --------
    :data:`CP2K_KEYS_ALIAS` : :class:`dict` [:class:`str`, :class:`tuple` [:class:`str`, ...]]
        A dictionary mapping ``key_path`` aliases to the actual keys.

    """  # noqa: E501
    if isinstance(value, Sequence):
        for prm_map in value:  # type: MappingScalar | MappingSequence
            set_prm(settings, key, prm_map, mol)
        return
    else:
        prm_map = copy.copy(value)

    try:
        prm_key = cast("str | Sequence[str]", prm_map.pop('param'))
    except KeyError as ex:
        raise KeyError("'param' has not been specified") from ex
    else:  # Extract the key path
        key_path = _get_key_path(key)

    # Get the list of settings located at **key_path**
    settings_base = settings.get_nested(key_path)
    if not isinstance(settings_base, list):  # Ensure it's a list of Settings
        settings_base = [settings_base] if settings_base else []
        settings.set_nested(key_path, settings_base)

    # charge, dipole and quadrupole is the only ff parameter assigned to a single atom,
    # rather than a sequence of n atoms
    atom_key = 'atom' if key in {'charge', 'dipole', 'quadrupole'} else 'atoms'

    # Map each pre-existing atom(-pair) to a list index in **prm_base**
    atom_map = {item.get(atom_key, None): i for i, item in enumerate(settings_base)}

    args = atom_map, settings_base, atom_key
    if isinstance(prm_key, Iterable):
        set_prm_values(prm_key, prm_map, *args)
        return
    raise TypeError(f"'param' expected a string or a sequence of strings; "
                    f"observed type: {prm_key.__class__.__name__!r}")


@overload
def set_prm_values(prm_key: str, prm_map: MappingScalar,
                   atom_map: MutableMapping[None | str, int],
                   settings_base: list[Settings], atom_key: str) -> None:
    ...
@overload  # noqa: E302
def set_prm_values(prm_key: Sequence[str], prm_map: MappingSequence,
                   atom_map: MutableMapping[None | str, int],
                   settings_base: list[Settings], atom_key: str) -> None:
    ...


def set_prm_values(prm_key, prm_map, atom_map,
                   settings_base, atom_key) -> None:  # noqa: E302
    """Assign the actual values specified in :func:`set_prm`.

    Parameters
    ----------
    prm_key : :class:`str` or :class:`~collections.abc.Sequence` [:class:`str`]
        The name(s) of the to-be set CP2K key(s), *e.g.* ``"sigma"`` and/or ``"epsilon"``.
        If ``iterable=False`` then this value should be a string;
        a sequence of strings is expected otherwise.

    prm_map : :class:`~collections.abc.MutableMapping` [:class:`str`, :data:`~typing.Any`]
        A dictionary containing the to-be set values.
        If ``iterable=False`` then its values should be scalars;
        sequences are expected otherwise.

    atom_map : :class:`~collections.abc.MutableMapping` [:class:`str`, :class:`int`]
        A dictionary for keeping track of which *atom_key* blocks are present in *settings_base*.

    settings_base : :class:`list` [:class:`qmflows.Settings`]
        A list of Settings to-be updated by *prm_map*.

    atom_key : :class:`str`
        The name of the CP2K atom key.
        Its value should be either ``"atoms"`` or ``"atom"`` depending on the parameter type
        of interest.

    """
    try:  # Read and parse the unit
        unit = prm_map.pop('unit')
    except KeyError:
        unit_iter = repeat('{}')
    else:
        unit_iter = _parse_unit(unit)  # type: ignore

    # Construct a DataFrame of parameters
    df = _construct_df(prm_key, prm_map)
    _validate_unit(unit_iter, df.columns)

    # Assign new parameters to the list of settings
    for atoms, prm in df.iterrows():
        prm_new = {k: unit.format(v) for unit, k, v in zip(unit_iter, df.columns, prm)}
        try:
            i: int = atom_map[atoms]
        except KeyError:
            s = Settings(prm_new)
            s[atom_key] = atoms
            settings_base.append(s)
            atom_map[atoms] = len(settings_base)
        else:
            settings_base[i].update(prm_new)


@singledispatch
def _parse_unit(unit: Iterable[None | str]) -> list[str]:
    """Convert *unit* into a to-be formatted string.

    *unit* can be eiher ``None``/ a string or an iterable consisting
    of aforementioned objects.

    """
    return [(f'[{u}] {{}}' if u is not None else '{}') for u in unit]


@_parse_unit.register(str)
@_parse_unit.register(type(None))
def _(unit: None | str) -> list[str]:
    if unit is None:
        return ['{}']
    else:
        return [f'[{unit}] {{}}']


@overload
def _construct_df(columns: Sequence[str], prm_map: MappingSequence) -> pd.DataFrame:
    ...
@overload  # noqa: E302
def _construct_df(columns: str, prm_map: MappingScalar) -> pd.DataFrame:
    ...
def _construct_df(columns, prm_map) -> pd.DataFrame:  # noqa: E302
    """Convert *prm_map* into a :class:`pandas.DataFrame` of strings with *columns* as columns.

    The main purpose of the DataFrame construction is to catch any errors
    due to an incorrect shape of either the *columns* or the *prm_map* values.

    *columns* and the *prm_map* values should be either both scalars or sequences.
    Scalars and sequences **cannot** be freely mixed.

    See :func:`set_prm_values`.

    """
    try:
        if isinstance(prm_map, NDFrame):
            data = prm_map.values.T
        else:
            data = np.array([v for v in prm_map.values()], dtype=np.str_)

        if data.ndim == 1:
            data.shape = -1, 1
            columns_ = [columns]
        else:
            columns_ = columns
        return pd.DataFrame(data, columns=columns_, index=prm_map.keys())

    except TypeError as ex:
        msg = str(ex)
        if msg.startswith('Index(...) must be called with a collection of some kind'):
            _raise_df_exc_scalar(columns, prm_map, ex)
        raise ex

    except ValueError as ex:
        _raise_df_exc_seq(columns, prm_map, ex)


def _raise_df_exc_seq(columns: Sequence[str], prm_map: MappingSequence, ex: Exception) -> NoReturn:
    """Raise an exception for :func:`_construct_df` using a sequence-based input."""
    msg = str(ex)
    column_count = len(columns)
    elements = 'elements' if column_count > 1 else 'element'

    # There are 2 types of ValueErrors which can be raised by incorrect column/values length
    if msg == 'setting an array element with a sequence':
        pass
    elif msg.startswith('Shape of passed values is '):
        pass
    elif msg.startswith('setting an array element with a sequence'):
        pass
    else:
        raise ex  # an Exception was raised due to some other unforseen reason

    try:
        for k, v in prm_map.items():
            if len(v) != column_count:  # Will raise a TypeError if it's an integer
                break
        else:
            raise ex  # This line should technically never be reached, but just in case

    except TypeError as ex2:
        raise TypeError(f"{k!r} expected a sequence with {column_count} {elements}; "
                        f"observed type: {v.__class__.__name__!r}") from ex2

    if not (isinstance(v, Sequence) or hasattr(v, '__array__')):
        raise TypeError(f"{k!r} expected a sequence with {column_count} {elements}; "
                        f"observed type: {v.__class__.__name__!r}") from ex

    raise LengthError(f"{k!r} expected a sequence with {column_count} {elements}; "
                      f"observed length: {len(v)}") from ex


def _raise_df_exc_scalar(columns: str, prm_map: MappingScalar, ex: Exception) -> NoReturn:
    """Raise an exception for :func:`_construct_df` using a scalar-based input."""
    # One of the values in *prm_map* is an iterable while it shouldn't be
    for k, v in prm_map.items():
        if isinstance(v, Iterable) and not isinstance(v, str):
            break
    else:
        raise ex  # This line should technically never be reached, but just in case

    raise TypeError(f"{k!r} expected a scalar; observed type: "
                    f"{v.__class__.__name__!r}") from ex


def _validate_unit(unit_iter: Iterator[str] | Sequence[str], columns: Sequence[str]) -> None:
    """Check if *unit_str* and *columns* in :func:`set_prm_values` are of the same length."""
    column_count = len(columns)
    if isinstance(unit_iter, Iterator):
        return  # It's a itertools.repeat instance; this is fine

    elif len(unit_iter) != column_count:
        elements = 'elements' if column_count > 1 else 'element'
        ex = ValueError(f"'unit' expected a sequence with {column_count} {elements}; "
                        f"observed length: {len(unit_iter)}")
        raise LengthError(ex) from ex


def _get_key_path(key: str | tuple[str, ...]) -> tuple[str, ...]:
    """Extract a value from :data:`CP2K_KEYS_ALIAS` if *key* is not a :class:`tuple`; return *key* otherwise.

    Id *key* is a :class:`tuple` it can have one of the three following structures:

    * The first key is ``"input"`` followed by the key path.
    * The first two keys are ``"specific"`` and ``"cp2k"`` followed by the key path.
    * *key* is equivalent to the key path, lacking any prepended keys.

    """  # noqa: E501
    if isinstance(key, tuple):  # It's a tuple with the key path alias
        if key[0] == 'input':
            return ('specific', 'cp2k') + key[1:]
        elif key[0:2] != ('specific', 'cp2k'):
            return ('specific', 'cp2k') + key
        else:
            return key

    try:  # It's an alias for a pre-defined key path (probably)
        return CP2K_KEYS_ALIAS[key]
    except KeyError as ex:
        raise KeyError(f"{key!r} section: no alias available for {key!r}") from ex


def prm_to_df(settings: MutableMapping) -> None:
    """Traverse the passed CP2K settings and convert all forcefield parameter blocks into instances of :class:`pandas.DataFrame`.

    Parameter blocks are recognized by their keys in one of the following two ways:
    * The key is also present in :data:`CP2K_KEYS_ALIAS`.
    * The key is a :class:`tuple`.

    Performs an inplace update of *settings*.

    Examples
    --------
    .. code:: python

        >>> from qmflows import Settings
        >>> from qmflows.cp2k_utils import prm_to_df

        >>> s = Settings()
        >>> s.lennard_jones = {
        ...     'param': ('epsilon', 'sigma'),
        ...     'unit': ('kcalmol', 'angstrom'),
        ...     'Cs': (1, 1),
        ...     'Cd': (2, 2),
        ...     'O': (3, 3),
        ...     'H': (4, 4)
        ... }

        >>> prm_to_df(s)

        >>> print(type(s.lennard_jones))
        <class 'pandas.core.frame.DataFrame'>

        >>> print(s)  # doctest: +NORMALIZE_WHITESPACE
        lennard_jones:                param      unit   Cs   Cd    O    H
                           epsilon  epsilon   kcalmol  1.0  2.0  3.0  4.0
                           sigma      sigma  angstrom  1.0  2.0  3.0  4.0


    Parameters
    ----------
    settings : :class:`qmflows.Settings`
        The CP2K input Settings

    See Also
    --------
    :data:`CP2K_KEYS_ALIAS`
        A dictionary mapping ``key_path`` aliases to the actual keys.

    """  # noqa: E501
    key_alias = CP2K_KEYS_ALIAS.keys()
    items = ((k, v) for k, v in settings.items() if k in key_alias or isinstance(k, tuple))
    for k, prm_map in items:
        try:
            columns = prm_map['param']
        except KeyError as ex:
            raise KeyError(f"{k} section: 'param' has not been specified") from ex

        df = _construct_df(columns, prm_map).T

        # All columns are now of the 'object' dtype; convert them into floats where possible
        for key, series in df.items():
            try:
                df[key] = series.astype(np.float64)
            except ValueError:  # Nop, not possible
                pass
        settings[k] = df


def _cp2k_keys_alias(indent: str = 8 * ' ') -> str:
    """Create a :class:`str` representations of :data:`CP2K_KEYS_ALIAS`.

    Used for constructing the module-level docstring.

    """
    width = 4 + max(len(k) for k in CP2K_KEYS_ALIAS)
    _mid = ',\n'.join(f'{(repr(k)+":"):{width}}{v!r}' for k, v in CP2K_KEYS_ALIAS.items())

    top = f'{indent}>>> CP2K_KEYS_ALIAS: dict[str, tuple[str, ...]] = ' + '{\n'
    mid = textwrap.indent(_mid, f'{indent}...     ')
    bot = f'\n{indent}... ' + '}'
    return f'{top}{mid}{bot}'


__doc__ = __doc__.format(cp2k_keys_alias=_cp2k_keys_alias())
