"""

Index
-----
.. currentmodule:: qmflows.packages.cp2k_mm
.. autosummary::
    CP2KMM
    CP2K_KEYS_ALIAS

API
---
.. autoclass:: CP2KMM
.. autodata:: CP2K_KEYS_ALIAS
    :annotation: : dict[str, tuple[str, ...]]

    .. code:: python

        >>> from typing import Dict, Tuple

{cp2k_keys_alias}

"""

import os
import textwrap
from os.path import join
from typing import Optional, Union, Any, Dict, ClassVar, Mapping, Sequence, Tuple
from warnings import warn

from scm import plams

from qmflows.parsers.cp2KParser import parse_cp2k_warnings
from qmflows.settings import Settings
from qmflows.warnings_qmflows import cp2k_warnings
from qmflows.packages.packages import (package_properties, parse_output_warnings, WarnMap)
from qmflows.packages.cp2k_package import CP2K_Result, CP2K

__all__ = ['cp2k_mm']


class CP2KMM(CP2K):
    """This class setup the requirement to run a CP2K Job <https://www.cp2k.org/>.

    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    **cp2k** function should be called.

    """

    generic_dict_file: ClassVar[str] = 'generic2CP2K.yaml'

    def __init__(self) -> None:
        super().__init__()

    def prerun(self, job_settings: Settings, mol: plams.Molecule,
               symbol_map: Optional[Mapping[str, str]],
               **kwargs: Any) -> None:
        """Run a set of tasks before running the actual job."""
        if symbol_map is None:
            symbol_map = {at.symbol: at.symbol for at in mol}
        _set_kinds(job_settings, symbol_map)

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_name: str = 'cp2k_job',
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> 'CP2K_Result':
        """Call the Cp2K binary using plams interface.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the numerical results.
        :type hdf5_file: String
        :param input_file_name: Optional name for the input.
        :type input_file_name: String
        :param out_file_name: Optional name for the output.
        :type out_file_name: String
        :param store_in_hdf5: wether to store the output arrays in HDF5 format.
        :type store_in_hdf5: Bool

        """
        # Input modifications
        cp2k_settings = Settings()
        cp2k_settings.input = settings.specific.cp2k

        # Create a Plams job
        job = plams.Cp2kJob(name=job_name, settings=cp2k_settings, molecule=mol)
        r = job.run()

        work_dir = work_dir if work_dir is not None else job.path

        warnings = parse_output_warnings(job_name, r.job.path,
                                         parse_cp2k_warnings, cp2k_warnings)

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        return CP2K_Result(cp2k_settings, mol, job_name, r.job.path, dill_path,
                           work_dir=work_dir, status=job.status, warnings=warnings)

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Create the settings input for complex cp2k keys.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param key: Special key declared in ``settings``.
        :param value: Value store in ``settings``.
        :param mol: molecular Geometry
        :type mol: plams Molecule

        """
        funs = {'psf': _parse_psf,
                'prm': _parse_prm}
        for k in CP2K_KEYS_ALIAS.keys():
            funs[k] = _set_prm

        # Function that handles the special keyword
        try:
            f = funs[key]
        except KeyError:  # Plan B: fall back to the CP2K super-class
            super().handle_special_keywords(settings, key, value, mol)
        else:
            f(settings, value, mol, key)


def _parse_psf(settings: Settings, key: str,
               value: Any, mol: plams.Molecule) -> None:
    """Assign a .psf file."""
    settings.specific.cp2k.force_eval.subsys.topology.conn_file_format = 'PSF'
    settings.specific.cp2k.force_eval.subsys.topology.conn_file_name = value


def _parse_prm(settings: Settings, key: str,
               value: Any, mol: plams.Molecule) -> None:
    """Assign a CHARMM-style .prm file."""
    settings.specific.cp2k.force_eval.mm.forcefield.parmtype = 'CHM'
    settings.specific.cp2k.force_eval.mm.forcefield.parm_file_name = value


def _set_kinds(s: Settings, symbol_map: Mapping[str, str]) -> None:
    """Generate the kind section for cp2k."""
    subsys = s.specific.cp2k.force_eval.subsys
    for custom_symbol, symbol in symbol_map.items():
        subsys[f'kind {custom_symbol}'].element = symbol


def _set_prm(settings: Settings, key: str,
             value: Mapping[str, Union[str, Sequence[str], float]],
             mol: plams.Molecule) -> None:
    """Assign a set of forcefield to **settings**.

    Example value for **value**.

    .. code:: python

        >>> print(value)
        key_path: [specific, cp2k, force_eval, mm, forcefield, nonbonded, lennard-jones]  # Mandatory
        unit: 'kcalmol'  # Optional
        Cs: 1.0  # The actual parameters
        Cd: 1.5
        O: 0.6666
        H: 1

    """
    # Read and parse the unit
    try:
        unit = value.pop('unit')
    except KeyError:
        unit_str = '{}'
    else:
        unit_str = f'[{unit}] {{}}'

    # Read and parse the key path; ensure the first two items are "specific" and "cp2k"
    key_path = _parse_key_path(value, key)

    # Get the list of settings located at **key_path**
    prm_base = settings.get_nested(key_path)
    if not isinstance(prm_base, list):  # Ensure it's a list of Settings
        prm_base = [prm_base]
        settings.set_nested(key_path, prm_base)

    # charge, dipole and quadrupole is the only ff parameter assigned to a single atom,
    # rather than a sequence of n atoms
    atom_key = 'atom' if key in {'charge', 'dipole', 'quadrupole'} else 'atoms'

    # Map each pre-existing atom(-pair) to a list index in **prm_base**
    atom_map = {item.get(atom_key, None): i for i, item in enumerate(prm_base)}

    # Assign new parameters to the list of settings
    for atoms, prm in value.items():
        try:
            i = atom_map[atoms]
        except KeyError:
            prm_base.append(Settings(
                {atom_key: atoms, key: unit_str.format(prm)}
            ))
            atom_map[atoms] = len(prm_base)
        else:
            prm_base[i][key] = unit_str.format(prm)


_BASE_PATH = ('specific', 'cp2k', 'force_eval', 'mm', 'forcefield')

#: A dictionary mapping ``keys`` aliases to the actual keys.
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


def _parse_key_path(value: Mapping[str, Union[str, Sequence[str]]],
                    prm_name: str) -> Sequence[str]:
    try:
        key_tup = value.pop('keys')
    except KeyError as ex:
        raise RuntimeError(f"{prm_name!r} section: 'keys' has not been specified") from ex

    if isinstance(key_tup, str):
        try:
            return CP2K_KEYS_ALIAS[key_tup]
        except KeyError as ex:
            raise RuntimeError(f"{prm_name!r} section: no 'keys' alias available "
                               f"for {value['keys']!r}") from ex

    first_key = key_tup[0]
    ret = ['specific', 'cp2k']

    if first_key == 'specific':
        ret = key_tup
    elif first_key == 'input':  # Plams style input
        ret += key_tup[1:]
    else:
        ret += key_tup
    return ret


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


__doc__ = __doc__.format(cp2k_keys_alias=_cp2k_keys_alias())
cp2k_mm = CP2KMM()
