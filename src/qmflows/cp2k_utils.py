"""A module with utility function for CP2K.

Index
-----
.. currentmodule:: qmflows.cp2k_utils
.. autosummary::
    set_prm
    map_psf_atoms
    CP2K_KEYS_ALIAS

API
---
.. autofunction:: set_prm
.. autofunction:: map_psf_atoms
.. autodata:: CP2K_KEYS_ALIAS
    :annotation: : dict[str, tuple[str, ...]]

    .. code:: python

        >>> from typing import Dict, Tuple

{cp2k_keys_alias}

"""

import textwrap
from os import PathLike
from io import TextIOBase
from functools import singledispatch
from itertools import repeat, islice
from collections import abc
from typing import (Union, Optional, List, Dict, Tuple, overload, MutableMapping,
                    Sequence, Any, AnyStr, Iterable)

import pandas as pd

from scm import plams
from qmflows.settings import Settings
from qmflows.backports import nullcontext

__all__ = ['set_prm', 'map_psf_atoms', 'CP2K_KEYS_ALIAS']

_BASE_PATH = ('specific', 'cp2k', 'force_eval', 'mm', 'forcefield')

#: A dictionary mapping ``key_path`` aliases to the actual keys.
CP2K_KEYS_ALIAS: Dict[str, Tuple[str, ...]] = {
    'bond': _BASE_PATH + ('bond',),
    'bend': _BASE_PATH + ('bend',),
    'ub': _BASE_PATH + ('bend', 'ub'),
    'torsion': _BASE_PATH + ('torsion',),
    'improper': _BASE_PATH + ('improper',),
    'charge': _BASE_PATH + ('charge',),
    'dipole': _BASE_PATH + ('dipole',),
    'quadrupole': _BASE_PATH + ('quadrupole',),

    'lennard-jones': _BASE_PATH + ('nonbonded', 'lennard-jones'),
    'lennard_jones': _BASE_PATH + ('nonbonded', 'lennard-jones'),
    'bmhft': _BASE_PATH + ('nonbonded', 'bmhft'),
    'bmhftd': _BASE_PATH + ('nonbonded', 'bmhftd'),
    'buck4ranges': _BASE_PATH + ('nonbonded', 'buck4ranges'),
    'buckmorse': _BASE_PATH + ('nonbonded', 'buckmorse'),
    'eam': _BASE_PATH + ('nonbonded', 'eam'),
    'genpot': _BASE_PATH + ('nonbonded', 'genpot'),
    'goodwin': _BASE_PATH + ('nonbonded', 'goodwin'),
    'ipbv': _BASE_PATH + ('nonbonded', 'ipbv'),
    'quip': _BASE_PATH + ('nonbonded', 'quip'),
    'siepmann': _BASE_PATH + ('nonbonded', 'siepmann'),
    'tersoff': _BASE_PATH + ('nonbonded', 'tersoff'),
    'williams': _BASE_PATH + ('nonbonded', 'williams'),

    'lennard-jones14': _BASE_PATH + ('nonbonded14', 'lennard-jones'),
    'lennard_jones14': _BASE_PATH + ('nonbonded14', 'lennard-jones'),
    'genpot14': _BASE_PATH + ('nonbonded14', 'genpot'),
    'goodwin14': _BASE_PATH + ('nonbonded14', 'goodwin'),
    'williams14': _BASE_PATH + ('nonbonded14', 'williams'),
}
del _BASE_PATH


MappingScalar = MutableMapping[str, Union[Optional[str], float]]
MappingSequence = MutableMapping[str, Union[Sequence[Optional[str]], Sequence[float]]]


class LengthError(ValueError):
    """A :exc:`ValueError` subclass for exceptions caused by incorrect lengths of a :class:`Mapping<collections.abc.Mapping>` or :class:`Sequence<collections.abc.Sequence>`."""  # noqa


def set_prm(settings: Settings, key: Union[str, Tuple[str, ...]],
            value: Union[MappingScalar, Sequence[MappingScalar], MappingSequence],
            mol: plams.Molecule) -> None:
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
    settings : :class:`qmflows.Settings<qmflows.settings.Settings>`
        The input CP2K settings.

    key : :class:`str` or :class:`tuple` [:class:`str`, ...]
        A path of CP2K keys or an alias for a pre-defined path (see :data:`CP2K_KEYS_ALIAS`).

    value : :class:`MutableMapping<collections.abc.MutableMapping>` [:class:`str`, ``T`` or :class:`Sequence<collections.abc.Sequence>` [``T``]]
        A dictionary with the to-be added parameters.
        Scalars and sequences **cannot** be freely mixed in the dictionary values;
        it should be one or the other.

        See the Examples section for more details.

    mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`
        A dummy argument in order to ensure signature compatiblity.

    Raises
    ------
    :exc:`RuntimeError`
        Raised when issues are encountered parsing either *key* or *value*.

    See Also
    --------
    :data:`CP2K_KEYS_ALIAS` : :class:`dict` [:class:`str`, :class:`tuple` [:class:`str`, ...]]
        A dictionary mapping ``key_path`` aliases to the actual keys.

    """  # noqa
    if isinstance(value, abc.Sequence):
        for prm_map in value:
            set_prm(settings, key, value, mol)
        return
    else:
        prm_map = value.copy()

    try:
        prm_key = prm_map.pop('param')
    except KeyError as ex:
        raise RuntimeError(f"{key!r} section: 'param' has not been specified") from ex
    else:  # Extract the key path
        key_path = _get_key_path(key)

    # Get the list of settings located at **key_path**
    settings_base = settings.get_nested(key_path)
    if not isinstance(settings_base, list):  # Ensure it's a list of Settings
        settings_base = [settings_base]
        settings.set_nested(key_path, settings_base)

    # charge, dipole and quadrupole is the only ff parameter assigned to a single atom,
    # rather than a sequence of n atoms
    atom_key = 'atom' if key in {'charge', 'dipole', 'quadrupole'} else 'atoms'

    # Map each pre-existing atom(-pair) to a list index in **prm_base**
    atom_map = {item.get(atom_key, None): i for i, item in enumerate(settings_base)}

    args = atom_map, settings_base, atom_key
    try:
        if isinstance(prm_key, abc.Iterable):
            return set_prm_values(prm_key, prm_map, *args)
        raise TypeError(f"'param' expected a string or a sequence of strings;\n"
                        f"observed type: {prm_key.__class__.__name__!r}")

    except LengthError as ex:
        raise RuntimeError(f'{key!r} section: {ex}') from ex.__cause__
    except Exception as ex:
        raise RuntimeError(f'{key!r} section: {ex}') from ex


@overload
def set_prm_values(prm_key: str, prm_map: MappingScalar,
                   atom_map: MutableMapping[Optional[str], int],
                   settings_base: List[Settings], atom_key: str) -> None:
    ...


@overload
def set_prm_values(prm_key: Sequence[str], prm_map: MappingSequence,
                   atom_map: MutableMapping[Optional[str], int],
                   settings_base: List[Settings], atom_key: str) -> None:
    ...


def set_prm_values(prm_key, prm_map, atom_map,
                   settings_base, atom_key) -> None:
    """Assign the actual values specified in :func:`set_prm`.

    Parameters
    ----------
    prm_key : :class:`str` or :class:`Sequence<collections.abc.Sequence>` [:class:`str`]
        The name(s) of the to-be set CP2K key(s), *e.g.* ``"sigma"`` and/or ``"epsilon"``.
        If ``iterable=False`` then this value should be a string;
        a sequence of strings is expected otherwise.

    prm_map : :class:`MutableMapping<collections.abc.MutableMapping>` [:class:`str`, ...]
        A dictionary containing the to-be set values.
        If ``iterable=False`` then its values should be scalars;
        sequences are expected otherwise.

    atom_map : :class:`MutableMapping<collections.abc.MutableMapping>` [:class:`str`, :class:`int`]
        A dictionary for keeping track of which *atom_key* blocks are present in *settings_base*.

    settings_base : :class:`list` [:class:`qmflows.Settings<qmflows.settings.Settings>`]
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
        unit_iter = _parse_unit(unit)

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
def _parse_unit(unit: Iterable[Optional[str]]) -> List[str]:
    """Convert *unit* into a to-be formatted string.

    *unit* can be eiher ``None``/ a string or an iterable consisting
    of aforementioned objects.

    """
    return [(f'[{u}] {{}}' if u is not None else '{}') for u in unit]


@_parse_unit.register(str)
@_parse_unit.register(type(None))
def _(unit: Optional[str]) -> List[str]:
    if unit is None:
        return ['{}']
    else:
        return [f'[{unit}] {{}}']


@singledispatch
def _construct_df(columns: Sequence[str], prm_map: MappingSequence) -> pd.DataFrame:
    """Convert *prm_map* into a :class:`pandas.DataFrame` of strings with *columns* as columns.

    The main purpose of the DataFrame construction is to catch any errors
    due to an incorrect shape of either the *columns* or the *prm_map* values.

    *columns* and the *prm_map* values should be either both scalars or sequences.
    Scalars and sequences **cannot** be freely mixed.

    See :func:`set_prm_values`.

    """
    try:
        return pd.DataFrame(prm_map, index=columns, dtype=str).T

    except ValueError as ex:
        msg = str(ex)
        column_count = len(columns)

        # There are 2 types of ValueErrors which can be raised by incorrect column/values length
        if msg == 'setting an array element with a sequence':
            pass
        elif msg.startswith(f'{column_count} columns passed, passed data had'):
            pass
        else:
            raise ex  # a ValueError was raised due to some other unforseen reason

        # One of the values in *prm_map* is of incorrect length; figure out which one
        for k, v in prm_map.items():
            if len(v) != column_count:
                break
        else:
            raise ex  # This line should technically never be reached, but just in case

        elements = 'elements' if column_count > 1 else 'element'
        raise LengthError(f"{k!r} expected a sequence with {column_count} {elements}; "
                          f"observed length: {len(v)}") from ex


@_construct_df.register(str)
def _(columns: str, prm_map: MappingScalar) -> pd.DataFrame:
    try:
        return pd.Series(prm_map, name=columns, dtype=str).to_frame()

    except ValueError as ex:
        msg = str(ex)

        # There is 1 type of ValueErrors which can be raised by incorrect column/values length
        if msg == 'setting an array element with a sequence':
            pass
        else:
            raise ex  # a ValueError was raised due to some other unforseen reason

        # One of the values in *prm_map* is an tierable while it shouldn't be
        for k, v in prm_map.items():
            if isinstance(v, abc.Iterable) and not isinstance(v, str):
                break
        else:
            raise ex  # This line should technically never be reached, but just in case

        raise TypeError(f"{k!r} expected a scalar; observed type: "
                        f"{v.__class__.__name__!r}") from ex


def _validate_unit(unit_iter: Union[repeat, Sequence[str]], columns: Sequence[str]) -> None:
    """Check if *unit_str* and *columns* in :func:`set_prm_values` are of the same length."""
    column_count = len(columns)
    if isinstance(unit_iter, repeat):
        return  # It's in infinite Iterator; this is fine

    elif len(unit_iter) != column_count:
        elements = 'elements' if column_count > 1 else 'element'
        ex = ValueError(f"'unit' expected a sequence with {column_count} {elements}; "
                        f"observed length: {len(unit_iter)}")
        raise LengthError(str(ex)) from ex


def _get_key_path(key: Union[str, Tuple[str, ...]]) -> Tuple[str, ...]:
    """Extract a value from :data:`CP2K_KEYS_ALIAS` if *key* is not a :class:`tuple`; return *key* otherwise.

    Id *key* is a :class:`tuple` it can have one of the three following structures:

    * The first key is ``"input"`` followed by the key path.
    * The first two keys are ``"specific"`` and ``"cp2k"`` followed by the key path.
    * *key* is equivalent to the key path, lacking any prepended keys.

    """  # noqa
    if isinstance(key, tuple):  # It's a tuple with the key path alias
        if key[0] == 'input':
            return ('specific', 'cp2k') + key[1:]
        elif not key[0:2] != ('specific', 'cp2k'):
            return ('specific', 'cp2k') + key
        else:
            return key

    try:  # It's an alias for a pre-defined key path (probably)
        return CP2K_KEYS_ALIAS[key]
    except KeyError as ex:
        raise RuntimeError(f"{key!r} section: no alias available for {key!r}") from ex


def map_psf_atoms(file: Union[AnyStr, PathLike, TextIOBase],
                  **kwargs: Any) -> Dict[str, str]:
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
    file : :class:`str`, :class:`PathLike<os.PathLike>` or :class:`TextIOBase<io.TextIOBase>`
        A path-like or file-like object containing the .psf file.
        Note that passed file-like objects should return strings (not bytes) upon iteration.

    /**kwargs : :data:`Any<typing.Any>`
        Further keyword arguments for :func:`open`.
        Only relevant when *file* is a path-like object.

    Returns
    -------
    :class:`dict` [:class:`str`, :class:`str`]
        A dictionary mapping atom types to atom names.
        Atom types/names are extracted from the passed .psf file.

    """  # noqa
    try:
        context_manager = open(file, **kwargs)  # path-like object
    except TypeError as ex:
        cls_name = file.__class__.__name__
        if cls_name == 'PSFContainer':  # i.e. the FOX.PSFContainer class
            return {k: v for v, k in zip(file.atom_type, file.atom_name)}
        elif not isinstance(file, abc.Iterator):
            raise TypeError("'file' expected a file- or path-like object; "
                            f"observed type: {cls_name!r}") from ex
        context_manager = nullcontext(file)  # a file-like object (hopefully)

    with context_manager as f:
        # A quick check to see *f* is not opened in bytes mode or something similar
        i = next(f)
        if not isinstance(i, str):
            raise TypeError(f"Iteration through {f!r} should yield a string; "
                            f"observed type: {i.__class__.__name__!r}")

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


def _cp2k_keys_alias(indent: str = f"{8 * ' '}") -> str:
    """Create a :class:`str` representations of :data:`CP2K_KEYS_ALIAS`.

    Used for constructing the module-level docstring.

    """
    width = 4 + max(len(k) for k in CP2K_KEYS_ALIAS)
    _mid = ',\n'.join(f'{(repr(k)+":"):{width}}{v!r}' for k, v in CP2K_KEYS_ALIAS.items())

    top = '{indent}>>> CP2K_KEYS_ALIAS: Dict[str, Tuple[str, ...]] = {\n'
    mid = textwrap.indent(_mid, f'{indent}...     ')
    bot = f'\n{indent}... ' + '}'
    return f'{top}{mid}{bot}'


__doc__ = __doc__.format(cp2k_keys_alias=_cp2k_keys_alias())
