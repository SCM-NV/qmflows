"""Deprecated alias for :mod:`qmflows.packages`."""

import warnings

from . import _settings as module
from .warnings_qmflows import QMFlowsDeprecationWarning

warnings.warn(
    f"`{__name__}` is a deprecated alias for `qmflows`",
    QMFlowsDeprecationWarning, stacklevel=2,
)

__globals__ = globals()
for k, v in vars(module).items():
    if k not in __globals__:
        __globals__[k] = v
del __globals__, k, v, module
