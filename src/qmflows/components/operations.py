"""Module containing some schedule utilities."""

__all__ = ['find_first_job', 'select_max', 'select_min']

from noodles import find_first, schedule
from scm.plams import Molecule

from qmflows.packages import Result
from qmflows.settings import Settings


@schedule
def find_first_job(
        pred: callable, packagelist: list, settings: Settings, molecule: Molecule,
        job_name: str, **kwargs) -> Result:
    """Return the first job to finish."""
    joblist = [package(
        settings, molecule, job_name=f"{package.pkg_name}_{job_name}", **kwargs)
        for package in packagelist]
    return find_first(pred, joblist)


@schedule
def select_max(results: Result, prop: str = 'energy'):
    """Select a result with the maximum value of a property from a Results list.

    If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    """
    filtered_results = [
        result for result in results if result and result.__getattr__(prop)]
    if len(filtered_results) > 0:
        selected = max(filtered_results, key=lambda item: getattr(item, prop))
        return selected
    else:
        return None


@schedule
def select_min(results: Result, prop='energy'):
    """Select a result with the minimum value of a property from a Results list.

    If the property is not available from a result (e.g. because
    the job failed) the result is ignored.
    """
    filtered_results = [
        result for result in results if result and result.__getattr__(prop)]
    if len(filtered_results) > 0:
        selected = min(filtered_results, key=lambda item: getattr(item, prop))
        return selected
    else:
        return None
