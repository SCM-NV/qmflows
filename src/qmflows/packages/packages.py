"""Common funcionality to call all the quantum packages."""

import fnmatch
import importlib
import inspect
import os
import sys
import warnings
from abc import abstractmethod, ABC
from types import ModuleType
from pathlib import Path
from functools import partial
from os.path import join
from warnings import warn
from typing import (
    Any, Callable, Optional, Union, ClassVar, Mapping, Iterator, Type, Dict, TypeVar, Tuple
)

import numpy as np
import pandas as pd
import pkg_resources as pkg
import scm.plams.interfaces.molecule.rdkit as molkit
from more_itertools import collapse
from noodles import has_scheduled_methods, schedule, serial
from noodles.run.threading.sqlite3 import run_parallel
from noodles.serial import AsDict, Registry
from noodles.serial.numpy import SerNumpyScalar, arrays_to_hdf5
from noodles.serial.path import SerPath
from noodles.serial.reasonable import SerReasonableObject
from rdkit import Chem
from scm import plams

from .serializer import SerMolecule, SerMol, SerSettings, SerNDFrame, SerReduce
from ..type_hints import WarnMap, WarnDict, WarnParser, PromisedObject, MolType, _Settings
from ..utils import InitRestart
from ..fileFunctions import yaml2Settings
from ..settings import Settings
from ..settings import _Settings as _SettingsType
from ..warnings_qmflows import QMFlows_Warning

__all__ = ['Package', 'run', 'registry', 'Result', 'SerMolecule', 'SerSettings']


def load_properties(name: str, prefix: str = 'properties') -> _Settings:
    """Load the properties-defining .yaml file from ."""
    file_name = os.path.join('data', 'dictionaries', f'{prefix}{name}.yaml')
    xs = pkg.resource_string("qmflows", file_name)
    return yaml2Settings(xs, mapping_type=_SettingsType)


class Result:
    """Class containing the results associated with a quantum chemistry simulation."""

    #: A :class:`Settings` instance with :class:`Result`-specific properties.
    #: Should be set when creating a subclass.
    prop_mapping: ClassVar[_Settings] = NotImplemented

    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 status: str = 'done',
                 warnings: Optional[WarnMap] = None) -> None:
        """Initialize a :class:`Result` instance.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param molecule: molecular Geometry
        :type molecule: plams Molecule
        :param job_name: Name of the computations
        :type job_name: str
        :param dill_path: The absolute path to the pickled .dill file.
        :type dill_path: str
        :param plams_dir: path to the ``Plams`` folder.
        :type plams_dir: str
        :param work_dir: scratch or another directory different from
        the `plams_dir`.
        type work_dir: str

        """
        plams_dir = None if plams_dir is None else Path(plams_dir)

        self.settings = settings
        self._molecule = molecule
        self.archive = {"plams_dir": plams_dir, 'work_dir': work_dir}
        self.job_name = job_name
        self.status = status
        self.warnings = warnings

        self._results_open = False
        self._results = dill_path

    def __deepcopy__(self, memo: Optional[Dict[int, Any]] = None) -> 'Result':
        """Return a deep copy of this instance."""
        cls = type(self)

        # Construct an empty instance while bypassing __init__()
        copy_instance = cls.__new__(cls)

        # Manually set all instance variables
        copy_instance.__dict__ = self.__dict__.copy()
        return copy_instance

    def __getattr__(self, prop: str) -> Any:
        """Return a section of the results.

        For example:

        ..code:: python

            >>> from qmflows.packages.packages import Result

            >>> result = Result(...)  # doctest: +SKIP
            >>> dipole = result.dipole  # doctest: +SKIP

        """
        is_private = prop.startswith('__') and prop.endswith('__')
        has_crashed = self.status in {'failed', 'crashed'}

        if not has_crashed and prop in self.prop_mapping:
            return self.get_property(prop)

        elif not (has_crashed or is_private or prop in self.prop_mapping):
            if self._results_open:
                warn(f"Generic property {prop!r} not defined",
                     category=QMFlows_Warning)

            # Do not issue this warning if the Results object is still pickled
            else:  # Unpickle the Results instance and try again
                self._unpack_results()
                return self.__getattr__(prop)

        elif has_crashed and not is_private:
            warn(f"""
            It is not possible to retrieve property: {prop!r}
            Because Job: {self.job_name!r} has {self.status}. Check the output.\n
            Are you sure that you have the package installed or
             you have loaded the package in the cluster. For example:
            `module load AwesomeQuantumPackage/3.141592`
            """, category=QMFlows_Warning)
        return None

    def get_property(self, prop: str) -> Any:
        """Look for the optional arguments to parse a property, which are stored in the properties dictionary."""  # noqa
        try:
            return super().__getattribute__(prop)
        except AttributeError:
            pass

        # Read the .yaml dictionary than contains the parsers names
        ds = self.prop_mapping[prop]

        # extension of the output file containing the property value
        file_ext = ds.get('file_ext')

        # If there is not work_dir returns None
        work_dir = self.archive.get('work_dir')

        # Plams dir
        plams_dir = self.archive['plams_dir']

        # Search for the specified output file in the folders
        if file_ext != "rkf":
            file_pattern = ds.get(
                'file_pattern', f'{self.job_name}*.{file_ext}')
        else:
            # AMS rename all the DFTB job names
            file_pattern = "dftb.rkf"

        output_files = list(collapse(map(partial(find_file_pattern, file_pattern),
                                         [plams_dir, work_dir])))
        if output_files:
            file_out = output_files[0]
            fun = getattr(import_parser(ds), ds['function'])

            # Read the keywords arguments from the properties dictionary
            kwargs = ds.get('kwargs', {})
            ret = ignore_unused_kwargs(fun, file_out, plams_dir=plams_dir, **kwargs)

            # Cache the property and return
            if sys.getsizeof(ret) < 10e5:
                setattr(self, prop, ret)
            return ret
        else:
            raise FileNotFoundError(f"""
            Property {prop} not found. No output file called: {file_pattern}. Folder used:
            plams_dir = {plams_dir}\n
            work_dir {work_dir}\n
            """)

    @property
    def results(self) -> Optional[plams.Results]:
        """Getter for :attr:`Result.results`.

        Get will load the .dill file and add all of its class attributes to this instance,
        barring the following three exceptions:

        * Private attributes/methods.
        * Magic methods.
        * Methods/attributes with names already contained within this instance.

        This attribute's value is set to ``None`` if the unpickling process fails.

        """
        if not self._results_open:
            self._unpack_results()
        return self._results

    def _unpack_results(self) -> None:
        """Unpack the pickled .dill file for :attr:`Results.results`."""
        self._results_open = True

        # Do not bother unpacking if None; i.e. if the job crashed
        if self._results is None:
            return

        # Ignore the Result.__getattr__() warnings for now
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=QMFlows_Warning)

            # Unpickle the results
            try:
                results = plams.load(self._results).results
                assert results is not None, f'Failed to unpickle {self._results!r}'
            except (AssertionError, plams.FileError) as ex:
                file_exc = ex
            else:
                file_exc = None
                attr_set = set(dir(self))
                for name in dir(results):
                    if name.startswith('_') or name in attr_set:
                        continue  # Skip methods which are either private, magic or preexisting

                    results_func = getattr(results, name)
                    setattr(self, name, results_func)

        # Failed to find or unpickle the .dill file; issue a warning
        if file_exc is not None:
            self._results = None
            warn(f"{file_exc}, setting value to 'None'",
                 category=QMFlows_Warning)
        else:
            self._results = results


PT = TypeVar("PT", bound="Package")


@has_scheduled_methods
class Package(ABC):
    """:class:`Package` is the base class to handle the invocation to different quantum package.

    The only relevant (instance) attribute of this class is :attr:`Package.pkg_name` which is a
    string representing the quantum package name that is going to be used to
    carry out the compuation.

    The life-cycle of :class:`Package` consists of 5 general steps:

    1. Initializing an instance: :meth:`Package.__init__`.
    2. Starting the job: :meth:`Package.__call__`.
       This method handles the task distribution between the instance's various methods.
    3. Converting all generic into specific settings: :meth:`Package.generic2specific`.
    4. Running the actual :class:`plams.Job<scm.plams.core.basejob.Job>`
       (including pre- and post-processing): :meth:`Package.run_job`.
    5. Returning the final :class:`Result` instance at the end of :meth:`Package.__call__`.

    """

    #: A class variable pointing to the :class:`Package`-specific :class:`Result` class.
    #: Should be set when creating a subclass.
    result_type: ClassVar[Type[Result]] = NotImplemented

    #: A class variable with the name of the generic .yaml file.
    #: Should be set when creating a subclass.
    generic_mapping: ClassVar[_Settings] = NotImplemented

    #: An instance variable with the name of the respective quantum chemical package.
    pkg_name: str

    def __init__(self, pkg_name: str) -> None:
        """Initialize a :class:`Package` instance.

        Parameters
        ----------
        pkg_name : :class:`str`
            The name of the respective quantum chemical package.
            See :attr:`Package.pkg_name`.

        """
        self.pkg_name = pkg_name

        # Ensure compatibility with the (typing-only) `builtins.function` class
        self.__name__: str = pkg_name
        self.__qualname__: str = pkg_name
        self.__annotations__: Dict[str, Any] = type(self).__call__.__annotations__

    def __reduce__(self: PT) -> Tuple[Type[PT], Tuple[str]]:
        """A helper function for :mod:`pickle`."""
        return type(self), (self.pkg_name,)

    @schedule(
        display="Running {self.pkg_name} {job_name}...",
        store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: MolType,
                 job_name: str = '', **kwargs: Any) -> Result:
        r"""Perform a job with the package specified by :attr:`Package.pkg_name`.

        Parameters
        ----------
        settings : :class:`~qmflows.settings.Settings`
            The user settings.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>` or :class:`rdkit.Mol<rdkit.Chem.rdchem.Mol>`
            A PLAMS or RDKit molecule to-be passed to the calculation.
        job_name : :class:`str`
            The name of the job.
        \**kwargs : :data:`~typing.Any`
            Further keyword arguments to-be passed to :meth:`Package.prerun`,
            :meth:`Package.run_job` and :meth:`Package.post_run`.

        Returns
        -------
        :class:`Result`
            A new Result instance.

        """  # noqa
        # Ensure that these variables have an actual value
        # Precaution against passing unbound variables to self.postrun()
        output_warnings = plams_mol = job_settings = None

        # There are not data from previous nodes in the dependecy trees
        # because of a failure upstream or the user provided None as argument
        if all(x is not None for x in [settings, mol]):
            #  Check if plams finishes normally
            try:
                # If molecule is an RDKIT molecule translate it to plams
                plams_mol = molkit.from_rdmol(
                    mol) if isinstance(mol, Chem.Mol) else mol

                if job_name != '':
                    kwargs['job_name'] = job_name

                # Settings transformations
                self.prerun(settings, plams_mol, **kwargs)
                job_settings = self.generic2specific(settings, mol)

                # Run the job
                result = self.run_job(job_settings, plams_mol, **kwargs)

                # Check if there are warnings in the output that render the calculation
                # useless from the point of view of the user
                warnings_tolerance = kwargs.get(
                    "terminate_job_in_case_of_warnings")
                output_warnings = result.warnings

                if None not in (warnings_tolerance, output_warnings):
                    issues = [w(msg) for msg, w in output_warnings.items()
                              if w in warnings_tolerance]
                    if issues:
                        warn(f"""
                        The Following Warning are rendered unacceptable in the Current
                        Workflow: {issues}\n
                        The results from Job: {job_name} are discarded.
                        """, category=QMFlows_Warning)
                        result = self.result_type(None, None, job_name=job_name,
                                                  dill_path=None, status='failed')

            # Otherwise pass an empty Result instance downstream
            except plams.core.errors.PlamsError as err:
                warn(f"Job {job_name} has failed.\n{err}",
                     category=QMFlows_Warning)
                result = self.result_type(None, None, job_name=job_name,
                                          dill_path=None, status='failed')
        else:
            warn(f"""
            Job {job_name} has failed. Either the Settings or Molecule
            objects are None, probably due to a previous calculation failure
            """, category=QMFlows_Warning)

            # Send an empty object downstream
            result = self.result_type(None, None, job_name=job_name,
                                      dill_path=None, status='failed')

        # Label this calculation as failed if there are not dependecies coming
        # from upstream
        self.postrun(result, output_warnings,
                     job_settings, plams_mol, **kwargs)
        return result

    def generic2specific(self, settings: Settings,
                         mol: Optional[plams.Molecule] = None) -> Settings:
        """Traverse *settings* and convert generic into package specific keys.

        Traverse all the key, value pairs of the *settings*, translating
        the generic keys into package specific keys as defined in the specific
        dictionary. If one key is not in the specific dictionary an error
        is raised. These new specific settings take preference over existing
        specific settings.

        Parameters
        ----------
        settings : :class:`~qmflows.settings.Settings`
            Settings provided by the user.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`, optional
            A PLAMS molecule to-be passed to the calculation.

        Returns
        -------
        :class:`~qmflows.settings.Settings`
            A new settings instance without any generic keys.

        """
        specific_from_generic_settings = Settings()
        for k, v in settings.items():
            if k == "specific":
                continue
            elif k == 'input':  # Allow for PLAMS-style input; i.e. settings.input.blablabla
                specific_from_generic_settings.specific[self.pkg_name].update(
                    v)
                continue

            if not self.generic_mapping.get(k):
                self.handle_special_keywords(
                    specific_from_generic_settings, k, v, mol)

        return settings.overlay(specific_from_generic_settings)

    def __repr__(self) -> str:
        """Create a string representation of this instance.

        Returns
        -------
        :class:`str`
            A string representation of this instnce.

        """
        vars_str = ', '.join(f'{k}={v!r}' for k,
                             v in sorted(vars(self).items()))
        return f'{self.__class__.__name__}({vars_str})'

    def prerun(self, settings: Settings, mol: plams.Molecule, **kwargs: Any) -> None:
        r"""Run a set of tasks before running the actual job.

        Parameters
        ----------
        settings : :class:`~qmflows.settings.Settings`
            Settings provided by the user.
            Note that these settings can still contain generic keywords.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`, optional
            A PLAMS molecule to-be passed to the calculation.
        \**kwargs : :data:`~typing.Any`
            Further keyword arguments to-be passed to :meth:`Package.run_job`.

        See Also
        --------
        :meth:`Package.run_job`
            A method which handles the running of
            the actual :class:`plams.Job<scm.plams.core.basejob.Job>`.

        """
        pass

    def postrun(self, result: Result,
                output_warnings: Optional[WarnMap] = None,
                settings: Optional[Settings] = None,
                mol: Optional[plams.Molecule] = None,
                **kwargs: Any) -> None:
        r"""Run a set of tasks after running the actual job.

        Parameters
        ----------
        result : :class:`Result`
            A Result instance.
        output_warnings : :class:`~collections.abc.Mapping` [:class:`str`, :class:`type` [:exc:`Warning`]], optional
            A Mapping which maps an error messages to Warning types.
        settings : :class:`~qmflows.settings.Settings`, optional
            User-provided Settings as processed by :meth:`Package.generic2specific`.
            Will be ``None`` if an error occured before this point.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`, optional
            A PLAMS molecule as passed to the calculation.
            Will be ``None`` if an error occured before
            the molecule was parsed in :meth:`Package.__call__`.
        \**kwargs : :data:`~typing.Any`
            Further keyword arguments that were passed to :meth:`Package.run_job`.

        See Also
        --------
        :meth:`Package.run_job`
            A method which handles the running of
            the actual :class:`plams.Job<scm.plams.core.basejob.Job>`.

        """  # noqa
        pass

    @staticmethod
    @abstractmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """`Abstract method <https://docs.python.org/3/library/abc.html#abc.abstractmethod>`_; should be implemented by the child class.

        A method providing additional processing for :class:`Package` dependant generic keywords.

        Parameters
        ----------
        settings : :class:`~qmflows.settings.Settings`, optional
            User-provided Settings as being processed by :meth:`Package.generic2specific`.
        key : :class:`str`
            The key associated with the special keyword
        value : :data:`~typing.Any`
            The value associated with the special *key*.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`
            A PLAMS molecule to-be passed to the calculation.

        See Also
        --------
        :meth:`Package.generic2specific`
            Traverse *settings* and convert generic into package specific keys.

        """  # noqa
        raise NotImplementedError("trying to call an abstract method")

    @classmethod
    @abstractmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule, job_name: str,
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> Result:
        r"""`Abstract method <https://docs.python.org/3/library/abc.html#abc.abstractmethod>`_; should be implemented by the child class.

        A method which handles the running of
        the actual :class:`plams.Job<scm.plams.core.basejob.Job>`.

        Parameters
        ----------
        settings : :class:`~qmflows.settings.Settings`, optional
            User-provided Settings as processed by :meth:`Package.generic2specific`.
        mol : :class:`plams.Molecule<scm.plams.mol.molecule.Molecule>`
            A PLAMS molecule to-be passed to the calculation.
        job_name : :class:`str`
            The name of the job.
        workdir : :class:`str` or :class:`~os.PathLike`, optional
            The path+folder name of the PLAMS working directory.
        \**kwargs : :data:`~typing.Any`
            Further keyword arguments.

        Returns
        -------
        :class:`Result`
            A new Result instance.

        """  # noqa
        raise NotImplementedError("The class representing a given quantum packages "
                                  "should implement this method")


def run(job: PromisedObject, runner: Optional[str] = None,
        path: Union[None, str, os.PathLike] = None,
        folder: Union[None, str, os.PathLike] = None,
        load_jobs: bool = False,
        **kwargs: Any) -> Result:
    r"""Pickup a runner and initialize it.

    Serves as a wrapper around :func:`noodles.run_parallel`.

    Parameters
    ----------
    job : :class:`noodles.PromisedObject<noodles.interface.PromisedObject>`
        The computation to run as constructed by :meth:`Package.__call__`.
    runner : :class:`str`, optional
        The job runner.
        Note that this value should be left at ``None``.
    path : :class:`str` or :class:`~os.PathLike`, optional
        The path where the PLAMS working directory will be created.
        Will default to the current working directory if ``None``.
    folder : :class:`str` or :class:`~os.PathLike`, optional
        The name of the new PLAMS working directory.
        Will default to ``"plams_workdir"`` if ``None``.
    load_jobs : :class:`bool`
        Load all pre-existing Jobs (contained within the working directory) into memory.
        Note that this can be quite slow if a large number of pre-existing jobs is present.
    \**kwargs : :data:`~typing.Any`
        Further keyword arguments to-be passed to :func:`call_default`.

    Returns
    -------
    :class:`Result`
        A new Result instance.
        The exact type will depend on **job**.

    See Also
    --------
    :func:`noodles.run_parallel`
        Run a workflow in parallel threads, storing results in a Sqlite3 database.

    """
    with InitRestart(path=path, folder=folder):
        plams.config.log.stdout = 0

        if runner is None:
            return call_default(job, kwargs.get('n_processes', 1),
                                kwargs.get('always_cache', True))
        else:
            raise ValueError(f"Don't know runner: {runner!r}")


def call_default(wf: PromisedObject, n_processes: int, always_cache: bool) -> Result:
    """Run locally using several threads.

    Caching can be turned off by specifying ``cache=None``.

    """
    # In case 'default_jobmanager' is not set (for some reason)
    try:
        workdir = plams.config.get('default_jobmanager').workdir
    except AttributeError as ex:
        raise plams.PlamsError("Failed to initialize the PLAMS jobmanager") from ex

    db_file = join(workdir, 'cache.db')
    return run_parallel(
        wf, n_threads=n_processes, registry=registry,
        db_file=db_file, always_cache=always_cache, echo_log=False)


#: A :class:`Registry` instance to-be returned by :func:`registry`.
REGISTRY: Registry = Registry(
    parent=serial.base() + arrays_to_hdf5(),
    types={
        Package: SerReduce(Package),
        Path: SerPath(),
        plams.Molecule: SerMolecule(),
        Chem.Mol: SerMol(),
        Result: AsDict(Result),
        Settings: SerSettings(),
        plams.KFFile: SerReasonableObject(plams.KFFile),
        plams.KFReader: SerReasonableObject(plams.KFReader),
        np.floating: SerNumpyScalar(),
        np.integer: SerNumpyScalar(),
        pd.DataFrame: SerNDFrame(pd.DataFrame),
        pd.Series: SerNDFrame(pd.Series),
    }
)


def registry() -> Registry:
    """Pass to the noodles infrastructure all the information related to the structure of the :class:`Package` object that is scheduled.

    This *Registry* class contains hints that help Noodles to encode
    and decode this Package object.

    """  # noqa
    return REGISTRY


def import_parser(ds: Mapping[str, str], module_root: str = "qmflows.parsers") -> ModuleType:
    """Import parser for the corresponding property."""
    module_sufix = ds['parser']
    module_name = module_root + '.' + module_sufix
    return importlib.import_module(module_name)


def find_file_pattern(path: Union[str, os.PathLike],
                      folder: Union[None, str, os.PathLike] = None) -> Iterator[str]:
    if folder is not None and os.path.exists(folder):
        return map(lambda x: join(folder, x), fnmatch.filter(os.listdir(folder), str(path)))
    else:
        return iter([])


def ignore_unused_kwargs(fun: Callable, *args: Any, **kwargs: Any) -> Any:
    """Inspect the signature of function `fun` and filter the keyword arguments.

    Searches for the keyword arguments which have a nonempty default values
    and then the dict `kwargs` those key-value pairs ignoring the rest.

    """
    ps = inspect.signature(fun).parameters

    # Look for the arguments with the nonempty defaults.
    defaults = filter(lambda t: t[1].default != inspect._empty, ps.items())

    # *kwargs* may contain keyword arguments not supported by *fun*
    # Extract from *kwargs* only the used keyword arguments
    kwargs2 = {k: kwargs[k] for k, _ in defaults if k in kwargs}
    return fun(*args, **kwargs2)


def parse_output_warnings(job_name: str,
                          plams_dir: Union[None, str, os.PathLike],
                          parser: WarnParser,
                          package_warnings: WarnMap) -> Optional[WarnDict]:
    """Look out for warnings in the output file."""
    output_files = find_file_pattern('*out', plams_dir)
    try:
        return parser(next(output_files), package_warnings)
    except StopIteration:
        warn(f"job: {job_name} has failed. check folder: {plams_dir}",
             category=QMFlows_Warning)
        return None
