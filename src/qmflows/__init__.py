"""QMFlows API."""

import sys
from typing import TYPE_CHECKING

from .__version__ import __version__

from .logger import logger

from .utils import InitRestart

from .packages import (
    adf, cp2k, cp2k_mm, dftb, orca, run, PackageWrapper)

from . import templates
from .settings import Settings

try:
    import rdkit
except ModuleNotFoundError as ex:
    _RDKIT_EX: "None | ModuleNotFoundError" = ex
else:
    _RDKIT_EX = None

__all__ = [
    '__version__',
    'logger',
    'InitRestart',
    'Angle', 'Dihedral', 'Distance', 'Settings',
    'adf', 'cp2k', 'cp2k_mm', 'dftb', 'orca', 'run', 'PackageWrapper',
    'example_H2O2_TS', 'example_freqs', 'example_generic_constraints',
    'example_partial_geometry_opt',
    'freq', 'geometry', 'singlepoint', 'ts', 'md', 'cell_opt',
    'find_first_job', 'select_max', 'select_min',
]

# Use `__getattr__` for loading (and copying) the templates in python >= 3.7
if TYPE_CHECKING or sys.version_info < (3, 7):
    from .templates import freq, geometry, singlepoint, ts, md, cell_opt

# Use `__getattr__` to raise a more descriptive error if RDKit
# is not installed (requires python >= 3.7)
if TYPE_CHECKING or sys.version_info < (3, 7) or _RDKIT_EX is None:
    from .components import (
        Angle,
        Dihedral,
        Distance,
        find_first_job,
        select_max,
        select_min,
    )
    from .examples import (
        example_H2O2_TS,
        example_freqs,
        example_generic_constraints,
        example_partial_geometry_opt,
    )
    from . import components, examples

if not TYPE_CHECKING and sys.version_info >= (3, 7):
    from ._init_utils import (
        getattr_method as __getattr__,
        dir_method as __dir__,
        RDKIT_SET,
    )
    if _RDKIT_EX is not None:
        __all__ = [name for name in __all__ if name not in RDKIT_SET]
    del RDKIT_SET
else:
    # Initalize the sub-module such that `_RDKIT_EX` can enter its namespace
    from . import _init_utils

# Clean up the namespace
del sys, TYPE_CHECKING, _RDKIT_EX
