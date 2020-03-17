"""A module with backports of objects added after Python 3.6."""

from contextlib import AbstractContextManager
from typing import Hashable

__all__ = ['nullcontext', 'Literal']

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

        def __init__(self, enter_result=None):
            self.enter_result = enter_result

        def __enter__(self):
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
            def __getitem__(self, name: Hashable) -> type:
                return type(name) if not isinstance(name, type) else name

        # Plan C; Literal.__getitem__ will now simply return the type
        # of the passed object; for example: Literal[True] == bool
        Literal = _Literal()
