"""A module with backports of objects added after Python 3.6."""

import sys
from contextlib import AbstractContextManager
from typing import Any, TypeVar, Optional, TYPE_CHECKING

__all__ = ['nullcontext']

# T, Final and Literal should be imported from qmflows.type_hints
T = TypeVar('T')


class _LiteralBackup:
    """A runtime-only placeholder for :class:`typing.Literal`."""

    def __getitem__(self, name):
        return Any


class _FinalBackup:
    """A runtime-only placeholder for :class:`typing.Final`."""

    def __getitem__(self, name):
        return Any


class _NullContextBackup(AbstractContextManager):
    """Context manager that does no additional processing.

    Used as a stand-in for a normal context manager, when a particular
    block of code is only sometimes used with a normal context manager:

    .. testsetup:: python

        >>> condition = False

    .. code:: python

        >>> cm = optional_cm if condition else nullcontext(1)
        >>> with cm as f:
        ...     print(f)
        1

    """

    def __init__(self, enter_result: T = None) -> None:
        self.enter_result = enter_result

    def __enter__(self) -> Optional[T]:
        return self.enter_result

    def __exit__(self, exc_type, exc_value, traceback):
        pass


# nullcontext was added in python 3.7
if sys.version_info >= (3, 7):
    from contextlib import nullcontext
else:
    nullcontext = _NullContextBackup
    nullcontext.__name__ = nullcontext.__qualname__ = 'nullcontext'


# Literal and Final were added to Python in 3.8;
# they were previously available in typing_extensions
if TYPE_CHECKING:
    if sys.version_info >= (3, 8):
        from typing import Literal, Final
    else:
        from typing_extensions import Literal, Final
else:
    _LiteralBackup.__name__ = _LiteralBackup.__qualname__ = 'Literal'
    _FinalBackup.__name__ = _FinalBackup.__qualname__ = 'Final'
    Literal = _LiteralBackup()
    Final = _FinalBackup()
