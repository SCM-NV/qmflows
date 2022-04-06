"""Deprecated alias for :mod:`qmflows.parsers`."""

import warnings

from . import _xyz as module

warnings.warn(
    f"`{__name__}` is a deprecated alias for `qmflows.parsers`",
    DeprecationWarning, stacklevel=2,
)

__globals__ = globals()
for k, v in vars(module).items():
    if k not in __globals__:
        __globals__[k] = v
del __globals__, k, v, module
