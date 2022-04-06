"""Deprecated alias for :mod:`qmflows.parsers.adf`."""

import warnings

from . import adf as module

warnings.warn(
    f"`{__name__}` is a deprecated alias for `{module.__name__}`",
    DeprecationWarning, stacklevel=2,
)

__globals__ = globals()
for k, v in vars(module).items():
    if k not in __globals__:
        __globals__[k] = v
del __globals__, k, v, module
