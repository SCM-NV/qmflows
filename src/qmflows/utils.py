"""Miscellaneous utilities."""

from __future__ import annotations

import os
import types
import shutil
import tempfile
from collections import Counter
from contextlib import AbstractContextManager, redirect_stdout
from functools import wraps
from os.path import abspath, join, normpath, splitext
from pathlib import Path
from contextlib import nullcontext
from collections.abc import Callable, Iterable, Iterator
from typing import IO, Any, TypeVar, TYPE_CHECKING

from scm.plams import JobManager, config, finish, init, load_all

from .type_hints import PathLike

if TYPE_CHECKING:
    _FT = TypeVar("_FT", bound=Callable[..., Any])

__all__ = ['to_runtime_error', 'init_restart', 'InitRestart']


def get_tmpfile_name(name: str = "tmp_qmflows_") -> Path:
    """Create a temporal file name."""
    return Path(tempfile.mkstemp(prefix=name)[1])


def to_runtime_error(func: _FT) -> _FT:
    """Decorate a `specific` function, translating any Exceptions into a :exc:`RuntimeError`.

    The Exception message is furthermore prepended with *key*.

    Examples
    --------
    .. code:: python

        >>> def func1(settings, key, value, mol):
        ...     raise Exception('error')

        >>> @to_runtime_error
        ... def func2(settings, key, value, mol):
        ...     raise Exception('error')

        >>> func1(None, 'func1', None, None)
        Traceback (most recent call last):
          ...
        Exception: error

        >>> func2(None, 'func2', None, None)
        Traceback (most recent call last):
          ...
        RuntimeError: 'func2' section: error

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
    return wrapper  # type: ignore[return-value]


def file_to_context(file: int | PathLike | IO[Any],
                    require_iterator: bool = True,
                    **kwargs: Any) -> AbstractContextManager[IO[Any]]:
    r"""Take a path- or file-like object and return an appropiate context manager instance.

    Passing a path-like object will supply it to :func:`open`,
    while passing a file-like object will pass it to :class:`contextlib.nullcontext`.

    Examples
    --------
    .. testsetup:: python

        >>> from pathlib import Path

        >>> path_like = Path('test') / 'test_files' / 'mypy.ini'

    .. code:: python

        >>> from io import StringIO

        >>> path_like = 'file_name.txt'  # doctest: +SKIP
        >>> file_like = StringIO('this is a file-like object')

        >>> context1 = file_to_context(path_like)
        >>> context2 = file_to_context(file_like)

        >>> with context1 as f1, context2 as f2:
        ...     pass # insert operations here


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
        return open(file, **kwargs)  # type: ignore

    # a file-like object (hopefully)
    except TypeError as ex:
        if require_iterator and not isinstance(file, Iterator):
            raise TypeError("'file' expected a file- or path-like object; "
                            f"observed type: {file.__class__.__name__!r}") from ex
        return nullcontext(file)  # type: ignore


def init_restart(
    path: None | str | os.PathLike[str] = None,
    folder: None | str | os.PathLike[str] = None,
    load_jobs: bool = False,
) -> None:
    """Call the PLAMS |init| function without creating a new directory.

    All pre-existing Jobs contained therein can be automatically loaded (see |load_all|)
    by setting the *load_jobs* keyword to ``True``.

    .. |init| replace:: :func:`plams.init()<scm.plams.core.functions.init>`
    .. |load_all| replace:: :func:`plams.load_all()<scm.plams.core.functions.load_all>`

    """
    is_init: bool = config.init
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
    # Remove the freshly created workdir and change to the previously created one.
    #
    # This branch is only relevant the first time `plams.init` is called (`config.init is False`)
    elif not is_init:
        shutil.rmtree(jobmanager.workdir)

    # This can happen if `init_restart` is launched multiple times
    elif not os.path.isdir(workdir):
        os.mkdir(workdir)

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
        >>> with InitRestart(path):  # doctest: +SKIP
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

    def __init__(
        self,
        path: None | str | os.PathLike[str] = None,
        folder: None | str | os.PathLike[str] = None,
        otherJM: None | Iterable[JobManager] = None,
        load_jobs: bool = False,
    ) -> None:
        """Initialize the context manager, assign the path, folder and jobmanagers."""
        self.path = path
        self.folder = folder
        self.otherJM = otherJM
        self.load_jobs = load_jobs

    def __enter__(self) -> None:
        """Enter the context manager, call :func:`init_restart`."""
        init_restart(path=self.path, folder=self.folder,
                     load_jobs=self.load_jobs)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: types.TracebackType | None,
    ) -> None:
        """Exit the context manager, call :func:`finish<scm.plams.core.functions.finish>`."""
        finish(self.otherJM)
