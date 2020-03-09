import textwrap
from typing import Union, Mapping, Sequence, Optional, List, Dict, Tuple

from scm import plams
from qmflows.settings import Settings

__all__ = ['set_prm', 'CP2K_KEYS_ALIAS']

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
    'genpot14': _BASE_PATH + ('nonbonded14', 'genpot'),
    'goodwin14': _BASE_PATH + ('nonbonded14', 'goodwin'),
    'williams14': _BASE_PATH + ('nonbonded14', 'williams'),
}
del _BASE_PATH


def set_prm(settings: Settings, key: str,
            value: Mapping[str, Union[str, Sequence[str], float, Sequence[float]]],
            mol: plams.Molecule) -> None:
    """Assign a set of forcefield to **settings**.

    Example value for **value**.

    .. code:: python

        >>> print(value)
        param: epsilon  # Mandatory
        unit: kcalmol  # Optional
        Cs: 1.0  # The actual parameters
        Cd: 1.5
        O: 0.6666
        H: 1

    """  # noqa
    # Get the list of settings located at **key_path**
    prm_base = settings.get_nested(key)
    if not isinstance(prm_base, list):  # Ensure it's a list of Settings
        prm_base = [prm_base]
        settings.set_nested(key, prm_base)

    # charge, dipole and quadrupole is the only ff parameter assigned to a single atom,
    # rather than a sequence of n atoms
    atom_key = 'atom' if key in {'charge', 'dipole', 'quadrupole'} else 'atoms'

    # Map each pre-existing atom(-pair) to a list index in **prm_base**
    atom_map = {item.get(atom_key, None): i for i, item in enumerate(prm_base)}


def _set_prm_scalar(prm_dict: Mapping[str, Union[str, float]],
                    atom_map: Mapping[Optional[str], int],
                    prm_base: List[Settings], key: str, atom_key: str):
    # Read and parse the unit
    try:
        unit = prm_dict.pop('unit')
    except KeyError:
        unit_str = '{}'
    else:
        unit_str = f'[{unit}] {{}}'

    # Assign new parameters to the list of settings
    for atoms, prm in prm_dict.items():
        try:
            i = atom_map[atoms]
        except KeyError:
            prm_base.append(Settings(
                {atom_key: atoms, key: unit_str.format(prm)}
            ))
            atom_map[atoms] = len(prm_base)
        else:
            prm_base[i][key] = unit_str.format(prm)


def _set_prm_sequence(prm_dict: Mapping[str, Union[Sequence[str], Sequence[float]]],
                      atom_map: Mapping[Optional[str], int],
                      prm_base: List[Settings], key: Sequence[str], atom_key: str):
    iterator = zip(prm_dict, key)
    for prm_dict_, key_ in iterator:
        _set_prm_scalar(prm_dict_, atom_map,)


def _cp2k_keys_alias(indent: str = f"{8 * ' '}... {4 * ' '}") -> str:
    """Create a :class:`str` representations of :data:`CP2K_KEYS_ALIAS`.

    Use for constructing the module-level docstring.

    """
    width = 4 + max(len(k) for k in CP2K_KEYS_ALIAS.keys())
    _mid = ',\n'.join(f'{(repr(k)+":"):{width}}{v!r}' for k, v in CP2K_KEYS_ALIAS.items())

    mid = textwrap.indent(_mid, indent)
    top = '        >>> CP2K_KEYS_ALIAS: Dict[str, Tuple[str, ...]] = {\n'
    bot = f'\n{indent[:-4]}' + '}'
    return f'{top}{mid}{bot}'