"""Deprecated alias for :mod:`qmflows.parsers.utils`."""

import warnings

from . import utils as module
from ..warnings_qmflows import QMFlowsDeprecationWarning

warnings.warn(
    f"`{__name__}` is a deprecated alias for `{module.__name__}`",
    QMFlowsDeprecationWarning, stacklevel=2,
)

__globals__ = globals()
for k, v in vars(module).items():
    if k not in __globals__:
        __globals__[k] = v
del __globals__, k, v, module
