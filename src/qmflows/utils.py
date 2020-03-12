"""Miscellaneous utilities."""

__author__ = "Felipe Zapata"

__all__ = ['dict2Setting', 'settings2Dict', 'zipWith', 'zipWith3', 'to_runtime_error',
           'init_restart', 'InitRestart']

import os
import shutil
from functools import wraps
from os.path import abspath, normpath, join, splitext
from typing import Union, Optional, Iterable, Callable
from contextlib import redirect_stdout, AbstractContextManager
from collections import Counter

from pymonad import curry
from scm.plams import config, init, finish, JobManager, load_all

from .settings import Settings


@curry
def zipWith(f, xs, ys):
    """zipWith generalises zip by zipping with the function given as the first argument"""
    return [f(*rs) for rs in zip(xs, ys)]


@curry
def zipWith3(f, xs, ys, zs):
    """The zipWith3 function takes a function which combines three elements,
    as well as three lists and returns a list of their point-wise combination.
    """
    return [f(*rs) for rs in zip(xs, ys, zs)]


def settings2Dict(s):
    """Transform a Settings object into a dict."""
    d = {}
    for k, v in s.items():
        if not isinstance(v, Settings):
            d[k] = v
        else:
            d[k] = settings2Dict(v)

    return d


def dict2Setting(d):
    """Transform recursively a dict into a Settings object."""
    r = Settings()
    for k, v in d.items():
        if isinstance(v, dict):
            r[k] = dict2Setting(v)
        else:
            r[k] = v

    return r


def to_runtime_error(func: Callable) -> Callable:
    """Decorate a `specific` function, translating any Exceptions into a :exc:`RuntimeError`.

    The Exception message is furthermore prepended with *key*.

    Examples
    --------
    .. code:: python

        >>> def func1(settings, key, value, mol):
        ...     raise Exception('error')

        >>> @to_runtime_error
        >>> def func2(settings, key, value, mol):
        ...     raise Exception('error')

        >>> func1(None, 'func1', None, None)
        Exception('error')

        >>> func1(None, 'func2', None, None)
        Exception('"func2" section: error')

    """
    @wraps(func)
    def wrapper(settings, key, value, mol, **kwargs):
        try:
            return func(settings, key, value, mol, *kwargs)
        except Exception as ex:
            if isinstance(ex, RuntimeError):
                raise ex
            raise RuntimeError(f"{key!r} section: {ex}").with_traceback(ex.__traceback__) from ex
    return wrapper


def init_restart(path: Union[None, str, os.PathLike] = None,
                 folder: Union[None, str, os.PathLike] = None,
                 load_jobs: bool = False) -> None:
    """Call the PLAMS |init| function without creating a new directory.

    All pre-existing Jobs contained therein can be automatically loaded (see |load_all|)
    by setting the *load_jobs* keyword to ``True``.

    .. |init| replace:: :func:`plams.init()<scm.plams.core.functions.init>`
    .. |load_all| replace:: :func:`plams.load_all()<scm.plams.core.functions.load_all>`

    """
    with open(os.devnull, 'w') as f, redirect_stdout(f):  # Temporary supress printing
        init(path, folder)

    # Parse variables
    path_ = os.getcwd() if path is None else abspath(path)
    folder_ = 'plams_workdir' if folder is None else normpath(folder)
    workdir = join(path_, folder_)

    # The JobManager instance
    jobmanager = config.default_jobmanager

    # There is no previously existing workdir; move along
    if jobmanager.workdir == workdir:
        return

    # There is an actual preexisting workdir;
    # Remove the freshly created workdir and change to the previously created one
    else:
        shutil.rmtree(jobmanager.workdir)

    # Update the files and folders in the default JobManager
    jobmanager.foldername = folder_
    jobmanager.workdir = workdir
    jobmanager.logfile = join(workdir, 'logfile')
    jobmanager.input = join(workdir, 'input')

    # Update JobManager.names
    folder_iterator = (splitext(f)[0] for f in os.listdir(workdir))
    jobmanager.names = dict(Counter(folder_iterator))

    # Load all previously pickled .dill files into the JobManager
    # NOTE: This can be quite slow if a large number of (large) jobs is stored therein
    if load_jobs:
        load_all(workdir, jobmanager)


class InitRestart(AbstractContextManager):
    """A context manager wrapper around :func:`init_restart`.

    Examples
    --------
    .. code:: python

        >>> path = "path/to/my/workdir"
        >>> with InitRestart(path):
        ...     ...  # Run any PLAMS Jobs here


    See Also
    --------
    :func:`plams.init()<scm.plams.core.functions.init>`:
        Initialize PLAMS environment. Create global ``config`` and the default
        :class:`JobManager<scm.plams.core.jobmanager.JobManager>`.

    :func:`plams.finish()<scm.plams.core.functions.finish>`:
        Wait for all threads to finish and clean the environment.

    :func:`plams.load_all()<scm.plams.core.functions.load_all>`:
        Load all jobs from *path*.

    """

    def __init__(self, path: Union[None, str, os.PathLike] = None,
                 folder: Union[None, str, os.PathLike] = None,
                 otherJM: Optional[Iterable[JobManager]] = None,
                 load_jobs: bool = False) -> None:
        """Initialize the context manager, assign the path, folder and jobmanagers."""
        self.path = path
        self.folder = folder
        self.otherJM = otherJM
        self.load_jobs = load_jobs

    def __enter__(self) -> None:
        """Enter the context manager, call :func:`init_restart`."""
        init_restart(path=self.path, folder=self.folder, load_jobs=self.load_jobs)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager, call :func:`finish<scm.plams.core.functions.finish>`."""
        finish(self.otherJM)
