"""A :class:`Package` subclass for classical CP2K calculations.

Index
-----
.. currentmodule:: qmflows.packages.cp2k_mm
.. autosummary::
    CP2KMM

API
---
.. autoclass:: CP2KMM

"""

import os
from os.path import join
from typing import Union, Any, ClassVar, Mapping, Dict, Callable

from scm import plams

from qmflows.cp2k_utils import set_prm, _map_psf_atoms, CP2K_KEYS_ALIAS
from qmflows.parsers.cp2KParser import parse_cp2k_warnings
from qmflows.settings import Settings
from qmflows.warnings_qmflows import cp2k_warnings
from qmflows.packages.packages import parse_output_warnings
from qmflows.packages.cp2k_package import CP2K_Result, CP2K

__all__ = ['cp2k_mm']


def _parse_psf(settings: Settings, key: str,
               value: Any, mol: plams.Molecule) -> None:
    """Assign a .psf file."""
    subsys = settings.specific.cp2k.force_eval.subsys
    if value is None:
        subsys.topology.use_element_as_kind = '.TRUE.'
        return

    symbol_map = _map_psf_atoms(None, key, value, None)
    for custom_symbol, symbol in symbol_map.items():
        subsys[f'kind {custom_symbol}'].element = symbol
    subsys.topology.conn_file_format = 'PSF'
    subsys.topology.conn_file_name = value


def _parse_prm(settings: Settings, key: str,
               value: Any, mol: plams.Molecule) -> None:
    """Assign a CHARMM-style .prm file."""
    settings.specific.cp2k.force_eval.mm.forcefield.parmtype = 'CHM'
    settings.specific.cp2k.force_eval.mm.forcefield.parm_file_name = value


SpecialFunc = Callable[[Settings, str, Any, plams.Molecule], None]

#: A :class:`dict` mapping special keywords to the appropiate function.
SPECIAL_FUNCS: Dict[str, SpecialFunc] = {'psf': _parse_psf, 'prm': _parse_prm}
for k in CP2K_KEYS_ALIAS:
    SPECIAL_FUNCS[k] = set_prm


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

    def prerun(self, settings: Settings, mol: plams.Molecule, **kwargs: Any) -> None:
        """Run a set of tasks before running the actual job."""
        cp2k = settings.specific.cp2k
        if not cp2k.get('psf'):
            cp2k.psf = None

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
        # Function that handles the special keyword
        if isinstance(key, tuple):
            set_prm(settings, value, mol, key)
        else:
            try:
                f = SPECIAL_FUNCS[key]
                f(settings, value, mol, key)
            except KeyError:  # Plan B: fall back to the CP2K super-class
                super().handle_special_keywords(settings, key, value, mol)


def _set_kinds(s: Settings, symbol_map: Mapping[str, str]) -> None:
    """Generate the kind section for cp2k."""
    subsys = s.specific.cp2k.force_eval.subsys
    for custom_symbol, symbol in symbol_map.items():
        subsys[f'kind {custom_symbol}'].element = symbol


cp2k_mm = CP2KMM()
