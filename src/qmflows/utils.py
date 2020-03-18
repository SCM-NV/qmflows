"""Miscellaneous utilities."""

__author__ = "Felipe Zapata"

__all__ = ['to_runtime_error', 'init_restart', 'InitRestart']

import os
import shutil
from functools import wraps
from os.path import abspath, normpath, join, splitext
from typing import Union, Optional, Iterable, Callable, Any, IO, ContextManager
from contextlib import redirect_stdout, AbstractContextManager
from collections import Counter, abc

from scm.plams import config, init, finish, JobManager, load_all

from .settings import Settings
from .backports import nullcontext
from .type_hints import PathLike


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

            exc = RuntimeError(f"{key!r} section: {ex}")
            raise exc.with_traceback(ex.__traceback__) from ex
    return wrapper


def file_to_context(file: Union[PathLike, IO],
                    require_iterator: bool = True,
                    **kwargs: Any) -> ContextManager[IO]:
    r"""Take a path- or file-like object and return an appropiate context manager instance.

    Passing a path-like object will supply it to :func:`open`,
    while passing a file-like object will pass it to :class:`contextlib.nullcontext`.

    Examples
    --------
    .. code:: python

        >>> from io import StringIO

        >>> path_like = 'file_name.txt'
        >>> file_like = StringIO('this is a file-like object')

        >>> context1 = file_to_context(path_like)
        >>> context2 = file_to_context(file_like)

        >>> with context1 as f1, with context2 as f2:
        ...     ... # insert operations here


    Parameters
    ----------
    file : :class:`str`, :class:`bytes`, :class:`os.PathLike` or :class:`io.TextIOBase`
        A `path- <https://docs.python.org/3/glossary.html#term-path-like-object>`_ or
        `file-like <https://docs.python.org/3/glossary.html#term-file-object>`_ object.

    require_iterator : :class:`bool`
        If``True``, loosen the constraints on what constitutes a file-like object
        and allow *file* to-be an :class:`~collections.abc.Iterator`.

    /**kwargs : :data:`~typing.Any`
        Further keyword arguments for :func:`open`.
        Only relevant if *file* is a path-like object.

    Returns
    -------
    :func:`open` or :class:`~contextlib.nullcontext`
        An initialized context manager.
        Entering the context manager will return a file-like object.

    """
    # path-like object
    try:
        return open(file, **kwargs)

    # a file-like object (hopefully)
    except TypeError as ex:
        if require_iterator and not isinstance(file, abc.Iterator):
            raise TypeError("'file' expected a file- or path-like object; "
                            f"observed type: {file.__class__.__name__!r}") from ex
        return nullcontext(file)


def init_restart(path: Optional[PathLike] = None,
                 folder: Optional[PathLike] = None,
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

    def __init__(self, path: Optional[PathLike] = None,
                 folder: Optional[PathLike] = None,
                 otherJM: Optional[Iterable[JobManager]] = None,
                 load_jobs: bool = False) -> None:
        """Initialize the context manager, assign the path, folder and jobmanagers."""
        self.path = path
        self.folder = folder
        self.otherJM = otherJM
        self.load_jobs = load_jobs

    def __enter__(self) -> None:
        """Enter the context manager, call :func:`init_restart`."""
        init_restart(path=self.path, folder=self.folder,
                     load_jobs=self.load_jobs)

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        """Exit the context manager, call :func:`finish<scm.plams.core.functions.finish>`."""
        finish(self.otherJM)
