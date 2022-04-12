"""Deprecated alias for :mod:`qmflows.packages._serializer`."""

import warnings

from . import _serializer as module
from ..warnings_qmflows import QMFlowsDeprecationWarning

warnings.warn(
    f"`{__name__}` is deprecated and will be removed in the future",
    QMFlowsDeprecationWarning, stacklevel=2,
)

__globals__ = globals()
for k, v in vars(module).items():
    if k not in __globals__:
        __globals__[k] = v
del __globals__, k, v, module
