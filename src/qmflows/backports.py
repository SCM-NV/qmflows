"""A module with backports."""

from typing import Union
from contextlib import AbstractContextManager

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


try:  # Literal was added to the typing module in Python 3.8
    from typing import Literal
except ImportError:
    try:  # It was available in the third-party `typing_extensions` module prior to 3.8
        from typing_extensions import Literal
    except ImportError:  # Plan C
        class _Literal:
            """Return ``typing.Union[name, name]`` when calling ``.__getitem__(name)``."""

            def __getitem__(self, name):
                return Union[name, name]

        Literal = _Literal()
