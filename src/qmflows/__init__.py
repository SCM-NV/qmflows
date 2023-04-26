"""QMFlows API."""

# flake8: noqa: E402

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

from ._version import __version__ as __version__
from ._version_info import version_info as version_info

from ._logger import logger

from .utils import InitRestart

from .packages import (
    adf, cp2k, cp2k_mm, dftb, orca, run, PackageWrapper)

from . import templates
from ._settings import Settings

try:
    import rdkit
except ModuleNotFoundError as ex:
    _RDKIT_EX: None | ModuleNotFoundError = ex
else:
    _RDKIT_EX = None

__all__ = [
    'logger',
    'InitRestart',
    'Angle', 'Dihedral', 'Distance', 'Settings',
    'adf', 'cp2k', 'cp2k_mm', 'dftb', 'orca', 'run', 'PackageWrapper',
    'freq', 'geometry', 'singlepoint', 'ts', 'md', 'cell_opt',
    'find_first_job', 'select_max', 'select_min',
]

if TYPE_CHECKING or _RDKIT_EX is None:
    from .components import (
        Angle,
        Dihedral,
        Distance,
        find_first_job,
        select_max,
        select_min,
    )
    from . import components, examples

if TYPE_CHECKING:
    from .templates import (
        freq,
        geometry,
        singlepoint,
        ts,
        md,
        cell_opt,
    )
else:
    # Use `__getattr__` to raise a more descriptive error if RDKit is not installed
    from ._init_utils import (
        getattr_method as __getattr__,
        dir_method as __dir__,
        RDKIT_SET,
    )
    if _RDKIT_EX is not None:
        __all__ = [name for name in __all__ if name not in RDKIT_SET]
    del RDKIT_SET

# Clean up the namespace
del sys, TYPE_CHECKING, _RDKIT_EX
