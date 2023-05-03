"""``__getattr__`` and ``__dir__`` implementations for the main QMFlows namespace."""

from __future__ import annotations

import types
from typing import Any

import qmflows

__all__ = ["dir_method", "getattr_method", "TEMPLATE_DICT", "RDKIT_SET", "RDKIT_EX"]

RDKIT_EX = qmflows._RDKIT_EX

# Map template names to the (cached) template
TEMPLATE_DICT = types.MappingProxyType({
    k: getattr(qmflows.templates, k) for k in qmflows.templates.__all__
})

# Map RDKit-requiring objects to their namespace
RDKIT_SET = frozenset({
    "components",
    "Angle",
    "Dihedral",
    "Distance",
    "find_first_job",
    "select_max",
    "select_min",
    "examples",
})


def __getattr__(self: types.ModuleType, name: str) -> Any:
    """Ensure that templates are always copied and the RDKit functions are loaded lazilly."""
    # Always return a copy of the template, as inplace operations will otherwise
    # modify the original template in the qmflows namespace
    try:
        return TEMPLATE_DICT[name].copy()
    except KeyError:
        pass

    if name in RDKIT_SET:
        raise ImportError(f"{name!r} requires the optional RDKit package") from RDKIT_EX
    raise AttributeError(f"module {self.__name__!r} has no attribute {name!r}")


def __dir__(self: types.ModuleType) -> list[str]:
    """Manually insert the qmflows templates and RDKit functions into :func:`dir`."""
    try:
        return self._DIR_CACHE.copy()
    except AttributeError:
        pass

    cache_set = set(object.__dir__(qmflows)) | TEMPLATE_DICT.keys()
    if RDKIT_EX is None:
        cache_set |= RDKIT_SET

    cache = sorted(cache_set)
    setattr(self, "_DIR_CACHE", cache)
    return cache.copy()


# Alias the functions under a different names such that they don't
# trigger module-level `getattr`/`dir` calls
_getattr_func = __getattr__
_dir_func = __dir__
del __getattr__
del __dir__

getattr_method = types.MethodType(_getattr_func, qmflows)
dir_method = types.MethodType(_dir_func, qmflows)
