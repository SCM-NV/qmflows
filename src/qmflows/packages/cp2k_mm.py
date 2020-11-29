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
from os.path import join, abspath
from typing import Union, Any, ClassVar, Dict, Type

from scm import plams

from .packages import Result, parse_output_warnings, load_properties
from .cp2k_package import CP2K
from ..cp2k_utils import set_prm, _map_psf_atoms, CP2K_KEYS_ALIAS
from ..parsers.cp2KParser import parse_cp2k_warnings
from ..settings import Settings
from ..warnings_qmflows import cp2k_warnings
from ..type_hints import Generic2Special, Final, _Settings

__all__ = ['cp2k_mm']


class CP2KMM_Result(Result):
    """A class providing access to CP2KMM result."""

    prop_mapping: ClassVar[_Settings] = load_properties('CP2KMM', prefix='properties')


class CP2KMM(CP2K):
    """A Package subclass for running `CP2K Jobs <https://www.cp2k.org/>`_ for classical forcefield calculations.

    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    :data:`cp2k_mm` function should be called.

    """  # noqa: E501

    generic_mapping: ClassVar[_Settings] = load_properties('CP2KMM', prefix='generic2')
    result_type: ClassVar[Type[Result]] = CP2KMM_Result

    def __init__(self, pkg_name: str = "cp2k") -> None:
        super().__init__(pkg_name)

    def prerun(self, settings: Settings, mol: plams.Molecule, **kwargs: Any) -> None:
        """Run a set of tasks before running the actual job."""
        psf = settings.get('psf')
        if not psf:
            settings.psf = None

        # Fix this at some point in the future
        """
        from pathlib import Path

        # Identify the number of pre-existing jobs
        jm = plams.config.default_jobmanager
        i = 1 + sum(jm.names.values())

        # Figure out the working direcyory
        workdir = kwargs.get('workdir')
        if workdir is None:
            workdir = Path(jm.path) / jm.foldername
        else:
            workdir = Path(workdir)

        # Set psf to None if not specified; write it if it's a FOX.PSFContainer instance
        psf = settings.get('psf')
        if not psf:
            settings.psf = None
        elif psf.__class__.__name__ == 'PSFContainer':
            psf_name = workdir / f"{kwargs.get('job_name', 'cp2k_job')}.{i}.psf"
            psf.write(psf_name)
            settings.psf = str(psf_name)

        # Write it if it's a FOX.PRMContainer instance
        prm = settings.get('prm')
        if prm.__class__.__name__ == 'PRMContainer':
            prm_name = workdir / f"{kwargs.get('job_name', 'cp2k_job')}.{i}.prm"
            prm.write(prm_name)
            settings.prm = str(prm_name)
        """

    @classmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule,
                job_name: str = 'cp2k_job',
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> CP2KMM_Result:
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

        return cls.result_type(cp2k_settings, mol, job_name, dill_path=dill_path,
                               plams_dir=r.job.path, work_dir=work_dir,
                               status=job.status, warnings=warnings)

    @classmethod
    def handle_special_keywords(cls, settings: Settings, key: str,
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
            set_prm(settings, key, value, mol)
        else:
            try:
                f = cls.SPECIAL_FUNCS[key]
            except KeyError:  # Plan B: fall back to the CP2K super-class
                super().handle_special_keywords(settings, key, value, mol)
            else:
                f(settings, key, value, mol)

    #: A :class:`dict` mapping special keywords to the appropiate function.
    SPECIAL_FUNCS: ClassVar[Dict[str, Generic2Special]]

    @staticmethod
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

        subsys.topology.conn_file_name = abspath(value)
        subsys.topology.conn_file_format = 'PSF'

    @staticmethod
    def _parse_prm(settings: Settings, key: str,
                   value: Any, mol: plams.Molecule) -> None:
        """Assign a CHARMM-style .prm file."""
        forcefield = settings.specific.cp2k.force_eval.mm.forcefield
        forcefield.parm_file_name = abspath(value)
        if not forcefield.parmtype:
            forcefield.parmtype = 'CHM'

    @staticmethod
    def _parse_periodic(s: Settings, key: str, value: Any, mol: plams.Molecule) -> None:
        """Set the keyword for periodic calculations."""
        force_eval = s.specific.cp2k.force_eval
        force_eval.subsys.cell.periodic = value
        force_eval.mm.poisson.periodic = value

        ewald = force_eval.mm.poisson.ewald
        if value.lower() == 'none':
            ewald.ewald_type = value
        elif not ewald.ewald_type:
            ewald.ewald_type = 'EWALD'


CP2KMM.SPECIAL_FUNCS = {
    'psf': CP2KMM._parse_psf,
    'prm': CP2KMM._parse_prm,
    'periodic': CP2KMM._parse_periodic
}

for k in CP2K_KEYS_ALIAS:
    CP2KMM.SPECIAL_FUNCS[k] = set_prm

#: An instance of :class:`CP2KMM`.
#: Only one instance of this class should exist at any given momemt;
#: *i.e.* this value is a singleton.
cp2k_mm: Final[CP2KMM] = CP2KMM()
