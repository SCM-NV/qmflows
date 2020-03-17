"""A module with backports of objects added after Python 3.6."""

from contextlib import AbstractContextManager
from typing import Any, Union, overload, Tuple, TypeVar, Type

__all__ = ['nullcontext']

# T and Literal should be imported from qmflows.type_hints
T = TypeVar('T')


try:
    from contextlib import nullcontext
except ImportError:  # nullcontext was added in python 3.7
    class nullcontext(AbstractContextManager):
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

        def __enter__(self) -> T:
            return self.enter_result

        def __exit__(self, exc_type, exc_value, traceback):
            pass


try:
    # Plan A: literal was added in Python 3.8
    from typing import Literal

except ImportError:
    try:
        # Plan B: literal was previously available in a third party package
        from typing_extensions import Literal

    except ImportError:
        class _Literal:
            @staticmethod
            def to_type(item: Any) -> type:
                """Return ``typing.Type[item]`` if *item* is a :class:`type` instance; return ``type(item)`` otherwhise."""  # noqa
                return Type[item] if isinstance(item, type) else type(item)

            @overload
            def __getitem__(self, name: Tuple[Union[Type[T], T], ...]) -> Union[Type[T]]: ...

            @overload
            def __getitem__(self, name: Union[Type[T], T]) -> Type[T]: ...

            def __getitem__(self, name):
                if isinstance(name, tuple):
                    type_tup = tuple(self.to_type(i) for i in name)
                    return Union[type_tup]
                else:
                    return self.to_type(name)

        # Plan C; Literal.__getitem__ will now simply return the type
        # of the passed object; for example: Literal[True] == bool
        Literal = _Literal()
