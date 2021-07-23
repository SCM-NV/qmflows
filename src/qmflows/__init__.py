"""QMFlows API."""

import sys
import types
import importlib as _importlib
from typing import TYPE_CHECKING, Any

from .__version__ import __version__

from .logger import logger

from .utils import InitRestart

from .packages import (
    adf, cp2k, cp2k_mm, dftb, orca, run, PackageWrapper)

from . import templates
from .settings import Settings

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
    from .components import (
        Angle, Dihedral, Distance, find_first_job, select_max, select_min
    )
    from .examples import (
        example_H2O2_TS, example_freqs, example_generic_constraints, example_partial_geometry_opt
    )
    from . import components, examples
else:
    _TEMPLATES = frozenset(templates.__all__)
    _REQUIRES_RDKIT = types.MappingProxyType({
        "components": "qmflows.components",
        "Angle": "qmflows.components",
        "Dihedral": "qmflows.components",
        "Distance": "qmflows.components",
        "find_first_job": "qmflows.components",
        "select_max": "qmflows.components",
        "select_min": "qmflows.components",
        "examples": "qmflows.examples",
        "example_H2O2_TS": "qmflows.examples",
        "example_freqs": "qmflows.examples",
        "example_generic_constraints": "qmflows.examples",
        "example_partial_geometry_opt": "qmflows.examples",
    })

    _DIR_CACHE: "None | list[str]" = None

    def __getattr__(name: str) -> Any:
        """Ensure that the qmflows templates are always copied before returning."""
        if name in _TEMPLATES:
            return getattr(templates, name).copy()

        # Lazily load (and cache) the content of `qmflows.examples` and `
        # qmflows.components` in order to avoid directly importing RDKit
        module_name = _REQUIRES_RDKIT.get(name)
        if module_name is not None:
            globals()[module_name] = module = _importlib.import_module(module_name)
            globals()[name] = ret = getattr(module, name, module)
            return ret
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    def __dir__() -> "list[str]":
        """Manually insert the qmflows templates into :func:`dir`."""
        global _DIR_CACHE
        if _DIR_CACHE is None:
            _DIR_CACHE = list(globals()) + templates.__all__ + list(_REQUIRES_RDKIT)
            _DIR_CACHE.sort()
        return _DIR_CACHE

# Clean up the namespace
del sys, types, TYPE_CHECKING, Any
