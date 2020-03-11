"""A module with backports of objects added after Python 3.6."""

from contextlib import AbstractContextManager

__all__ = ['nullcontext']

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
