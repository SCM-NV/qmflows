"""Module containing some schedule utilities."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, TYPE_CHECKING

from noodles import find_first, schedule
from scm.plams import Molecule

from ..packages import Result, Package
from .._settings import Settings
from ..type_hints import PromisedObject

if TYPE_CHECKING:
    #: A function which takes a PromisedObject as argument and returns a boolean.
    Predicate = Callable[[PromisedObject], bool]

__all__ = ['find_first_job', 'select_max', 'select_min']


@schedule
def find_first_job(pred: Predicate,
                   packagelist: Iterable[Package],
                   settings: Settings,
                   molecule: Molecule,
                   job_name: str, **kwargs: Any) -> None | Result:
    """Return the first job to finish."""
    joblist = [package(
        settings, molecule, job_name=f"{package.pkg_name}_{job_name}", **kwargs)
        for package in packagelist]
    return find_first(pred, joblist)


@schedule
def select_max(results: Iterable[Result], prop: str = 'energy') -> None | Result:
    """Select a result with the maximum value of a property from a Results list.

    If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    """
    generator = (result for result in results if getattr(result, prop, None))
    try:
        return max(generator, key=lambda item: getattr(item, prop))
    except ValueError:  # Raised if generator is empty
        return None


@schedule
def select_min(results: Iterable[Result], prop: str = 'energy') -> None | Result:
    """Select a result with the minimum value of a property from a Results list.

    If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    """
    generator = (result for result in results if getattr(result, prop, None))
    try:
        return min(generator, key=lambda item: getattr(item, prop))
    except ValueError:  # Raised if generator is empty
        return None
