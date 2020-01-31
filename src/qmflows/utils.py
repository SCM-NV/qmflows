"""Miscellaneous utilities."""

__author__ = "Felipe Zapata"

__all__ = ['dict2Setting', 'settings2Dict', 'zipWith', 'zipWith3', 'init_restart', 'InitRestart']

import os
import shutil
from typing import Union, AnyStr, Optional, Iterable
from contextlib import redirect_stdout, AbstractContextManager

from pymonad import curry
from scm.plams import config, init, finish, JobManager

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


def init_restart(path: Union[None, AnyStr, os.PathLike] = None,
                 folder: Union[None, AnyStr, os.PathLike] = None) -> None:
    """Call the PLAMS |init| function without creating a new directory.

    .. |init| replace:: :func:`init<scm.plams.core.functions.init>`

    """
    with open(os.devnull, 'w') as f, redirect_stdout(f):  # Temporary supress printing
        init(path, folder)
    shutil.rmtree(config.default_jobmanager.workdir)  # Remove the freshly created workdir

    # Parse variables
    path_ = os.getcwd() if path is None else os.path.abspath(path)
    folder_ = 'plams_workdir' if folder is None else os.path.normpath(folder)
    workdir = os.path.join(path_, folder_)

    # Update the files and folders in the default JobManager
    config.default_jobmanager.foldername = folder_
    config.default_jobmanager.workdir = workdir
    config.default_jobmanager.logfile = os.path.join(workdir, 'logfile')
    config.default_jobmanager.input = os.path.join(workdir, 'input')


class InitRestart(AbstractContextManager):
    """A context manager wrapper around :func:`init_restart`."""

    def __init__(self, path: Union[None, AnyStr, os.PathLike] = None,
                 folder: Union[None, AnyStr, os.PathLike] = None,
                 otherJM: Optional[Iterable[JobManager]] = None) -> None:
        """Initialize the context manager, assign the path, folder and jobmanagers."""
        self.path = path
        self.folder = folder
        self.otherJM = otherJM

    def __enter__(self) -> None:
        """Enter the context manager, call :func:`init_restart`."""
        init_restart(self.path, self.folder)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager, call :func:`finish<scm.plams.core.functions.finish>`."""
        finish(self.otherJM)
