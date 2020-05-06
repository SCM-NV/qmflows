"""A module with backports of objects added after Python 3.6."""

import sys
from contextlib import AbstractContextManager
from typing import Any, Union, overload, Tuple, TypeVar, Type, Optional, TYPE_CHECKING

__all__ = ['nullcontext']

# T, Final and Literal should be imported from qmflows.type_hints
T = TypeVar('T')


class _LiteralBackup:
    @staticmethod
    def _to_type(item: Any) -> type:
        """Return ``typing.Type[item]`` if *item* is a :class:`type` instance; return ``type(item)`` otherwhise."""  # noqa
        return Type[item] if isinstance(item, type) else type(item)

    @overload
    def __getitem__(self, name: Tuple[Union[Type[T], T], ...]) -> Union[Type[T]]: ...
    @overload
    def __getitem__(self, name: Union[Type[T], T]) -> Type[T]: ...
    def __getitem__(self, name):
        if isinstance(name, tuple):
            type_tup = tuple(self._to_type(i) for i in name)
            return Union[type_tup]
        else:
            return self._to_type(name)


class _FinalBackup:
    def __getitem__(self, name: type) -> type:
        if not isinstance(name, type):
            raise TypeError(f'{self.__class__.__name__} accepts only single '
                            f'type Got {name!r:.100}.')
        return name


class _NullContextBackup(AbstractContextManager):
    """Context manager that does no additional processing.

    Used as a stand-in for a normal context manager, when a particular
    block of code is only sometimes used with a normal context manager:

    .. code:: python

        >>> cm = optional_cm if condition else nullcontext()
        >>> with cm:
        ...     ...  # Perform operation, using optional_cm if condition is True

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
