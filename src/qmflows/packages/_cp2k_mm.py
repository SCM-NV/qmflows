"""A :class:`~qmflows.packages.Package` subclass for classical CP2K calculations."""

from __future__ import annotations

import os
from warnings import warn
from os.path import abspath
from typing import Any, ClassVar, TYPE_CHECKING, Final

import numpy as np
from scm import plams

from ._packages import load_properties
from ._cp2k import CP2K, CP2K_Result
from ..cp2k_utils import set_prm, _map_psf_atoms, CP2K_KEYS_ALIAS
from .._settings import Settings
from ..warnings_qmflows import Key_Warning
from ..type_hints import _Settings

__all__ = ['CP2KMM_Result', 'CP2KMM', 'cp2k_mm']


class CP2KMM_Result(CP2K_Result):
    """A class providing access to CP2KMM result."""

    prop_mapping: ClassVar[_Settings] = load_properties('CP2KMM', prefix='properties')


class CP2KMM(CP2K):
    """A Package subclass for running `CP2K Jobs <https://www.cp2k.org/>`_ for classical forcefield calculations.

    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    :class:`~qmflows.cp2k_mm` function should be called.

    """  # noqa: E501

    generic_mapping: ClassVar[_Settings] = load_properties('CP2KMM', prefix='generic2')
    result_type: ClassVar[type[CP2KMM_Result]] = CP2KMM_Result

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

    if TYPE_CHECKING:
        @classmethod
        def run_job(
            cls,
            settings: Settings,
            mol: plams.Molecule,
            job_name: str = 'cp2k_job',
            work_dir: None | str | os.PathLike[str] = ...,
            validate_output: bool = True,
            **kwargs: Any,
        ) -> CP2KMM_Result:
            """Call the ORCA binary using plams interface."""
            ...

    @classmethod
    def handle_special_keywords(
        cls, settings: Settings, key: str, value: Any, mol: plams.Molecule
    ) -> None:
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
            return set_prm(settings, key, value, mol)

        f = cls.SPECIAL_FUNCS.get(key)
        if f is None:
            f = CP2K.SPECIAL_FUNCS.get(key)

        if f is None:
            warn(f'Generic keyword {key!r} not implemented for package CP2K', category=Key_Warning)
        else:
            f(settings, key, value, mol)

    @staticmethod
    def _parse_psf(settings: Settings, key: str,
                   value: None, mol: plams.Molecule) -> None:
        """Assign a .psf file."""
        subsys = settings.specific.cp2k.force_eval.subsys
        if value is None:
            symbol_list = sorted({at.symbol for at in mol})
            for symbol in symbol_list:
                subsys[f'kind {symbol}'].element = symbol
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
            ewald.ewald_type = 'SPME'

    @staticmethod
    def _parse_gmax(s: Settings, key: str, value: Any, mol: plams.Molecule) -> None:
        """Set the keyword for the ``gmax`` parameter."""
        ar = np.asarray(value)
        if not issubclass(ar.dtype.type, np.str_):
            ar = ar.astype(np.int64, copy=False, casting='same_kind')

        ewald = s.specific.cp2k.force_eval.mm.poisson.ewald
        if ar.ndim == 0:
            ewald.gmax = str(ar.item())
        elif ar.ndim == 1:
            ewald.gmax = " ".join(str(i) for i in ar)
        else:
            raise RuntimeError(f"gmax:{value!r}\nformat not recognized")


CP2KMM.SPECIAL_FUNCS = {
    'psf': CP2KMM._parse_psf,
    'prm': CP2KMM._parse_prm,
    'periodic': CP2KMM._parse_periodic,
    'gmax': CP2KMM._parse_gmax,
}

for k in CP2K_KEYS_ALIAS:
    CP2KMM.SPECIAL_FUNCS[k] = set_prm

cp2k_mm: Final[CP2KMM] = CP2KMM()
