"""QMFlows API."""

import sys
from typing import TYPE_CHECKING

from .__version__ import __version__

from .logger import logger

from .utils import InitRestart

from .components import (
    Angle, Dihedral, Distance,
    find_first_job, select_max, select_min)

from .packages import (
    adf, cp2k, cp2k_mm, dftb, orca, run, PackageWrapper)

from . import templates
from .settings import Settings
from .examples import (example_H2O2_TS, example_freqs, example_generic_constraints,
                       example_partial_geometry_opt)

__all__ = [
    '__version__',
    'logger',
    'InitRestart',
    'Angle', 'Dihedral', 'Distance', 'Settings',
    'adf', 'cp2k', 'cp2k_mm', 'dftb', 'orca', 'run', 'PackageWrapper',
    'example_H2O2_TS', 'example_freqs', 'example_generic_constraints',
    'example_partial_geometry_opt',
    'freq', 'geometry', 'singlepoint', 'ts', 'md', 'cell_opt',
    'find_first_job', 'select_max', 'select_min']

if TYPE_CHECKING or sys.version_info < (3, 7):
    from .templates import freq, geometry, singlepoint, ts, md, cell_opt
else:
    _TEMPLATES = frozenset(templates.__all__)
    _DIR_CACHE: "None | list[str]" = None

    def __getattr__(name: str) -> Settings:
        """Ensure that the qmflows templates are always copied before returning."""
        if name in _TEMPLATES:
            return getattr(templates, name).copy()
        raise AttributeError(f"module {__name__} has no attribute {name}")

    def __dir__() -> "list[str]":
        """Manually insert the qmflows templates into :func:`dir`."""
        global _DIR_CACHE
        if _DIR_CACHE is None:
            _DIR_CACHE = sorted(list(globals()) + templates.__all__)
        return _DIR_CACHE

# Clean up the namespace
del sys, TYPE_CHECKING
