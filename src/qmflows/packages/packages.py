"""Common funcionality to call all the quantum packages."""

import base64
import fnmatch
import importlib
import inspect
import os
import uuid
import warnings
from abc import abstractmethod, ABC
from types import ModuleType
from pathlib import Path
from functools import partial
from os.path import join
from warnings import warn
from typing import (Any, Callable, Optional, Dict, Union, ClassVar,
                    Iterable, Mapping, Iterator)

import numpy as np
import pkg_resources as pkg
import scm.plams.interfaces.molecule.rdkit as molkit
from more_itertools import collapse
from noodles import has_scheduled_methods, schedule, serial
from noodles.interface import PromisedObject
from noodles.run.threading.sqlite3 import run_parallel
from noodles.serial import AsDict, Registry, Serialiser
from noodles.serial.numpy import SerNumpyScalar, arrays_to_hdf5
from noodles.serial.path import SerPath
from noodles.serial.reasonable import SerReasonableObject
from rdkit import Chem
from scm import plams

from ..fileFunctions import yaml2Settings
from ..settings import Settings
from ..type_hints import WarnMap, WarnDict, WarnParser
from ..warnings_qmflows import QMFlows_Warning

__all__ = ['package_properties',
           'Package', 'run', 'registry', 'Result',
           'SerMolecule', 'SerSettings']

_BASE_PATH = Path('data') / 'dictionaries'

#: A dictionary mapping package names to .yaml files.
package_properties: Dict[Optional[str], Path] = {
    None: _BASE_PATH / 'propertiesNone.yaml',
    'adf': _BASE_PATH / 'propertiesADF.yaml',
    'dftb': _BASE_PATH / 'propertiesDFTB.yaml',
    'cp2k': _BASE_PATH / 'propertiesCP2K.yaml',
    'cp2kmm': _BASE_PATH / 'propertiesCP2KMM.yaml',
    'gamess': _BASE_PATH / 'propertiesGAMESS.yaml',
    'orca': _BASE_PATH / 'propertiesORCA.yaml'
}
del _BASE_PATH


class Result:
    """Class containing the results associated with a quantum chemistry simulation."""

    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 properties: Union[None, str, os.PathLike] = None,
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
        :param properties: path to the `yaml` file containing the data to
                           load the parser on the fly.
        :type properties: str

        """
        plams_dir = None if plams_dir is None else Path(plams_dir)
        self.settings = settings
        self._molecule = molecule
        xs = pkg.resource_string("qmflows", str(properties))
        self.prop_dict = yaml2Settings(xs)
        self.archive = {"plams_dir": plams_dir,
                        'work_dir': work_dir}
        self.job_name = job_name
        self.status = status
        self.warnings = warnings

        self._results_open = False
        self._results = dill_path

    def __deepcopy__(self, memo: Any) -> 'Result':
        """Return a deep copy of this instance."""
        cls = type(self)
        if not self._results_open or self._results is None:
            dill_path = self._results
        else:
            job = self._results.job
            dill_path = join(job.path, f"{job.name}.dill")

        return cls(self.settings,
                   self._molecule,
                   self.job_name,
                   dill_path,
                   plams_dir=self.archive['plams_dir'],
                   work_dir=self.archive['work_dir'],
                   status=self.status,
                   warnings=self.warnings
                   )

    def __getattr__(self, prop: str) -> Any:
        """Return a section of the results.

        For example:

        ..code:: python

            >>> from qmflows.packages.packages import Results

            >>> result = Result(...)
            >>> dipole = result.dipole

        """
        is_private = prop.startswith('__') and prop.endswith('__')
        has_crashed = self.status in {'failed', 'crashed'}

        if not has_crashed and prop in self.prop_dict:
            return self.get_property(prop)

        elif not (has_crashed or is_private or prop in self.prop_dict):
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
        # Read the .yaml dictionary than contains the parsers names
        ds = self.prop_dict[prop]

        # extension of the output file containing the property value
        file_ext = ds['file_ext']

        # If there is not work_dir returns None
        work_dir = self.archive.get('work_dir')

        # Plams dir
        plams_dir = self.archive['plams_dir']

        # Search for the specified output file in the folders
        file_pattern = ds.get('file_pattern', f'{self.job_name}*.{file_ext}')

        output_files = list(collapse(map(partial(find_file_pattern, file_pattern),
                                         [plams_dir, work_dir])))
        if output_files:
            file_out = output_files[0]
            fun = getattr(import_parser(ds), ds['function'])
            # Read the keywords arguments from the properties dictionary
            kwargs = ds.get('kwargs', {})
            kwargs['plams_dir'] = plams_dir
            return ignored_unused_kwargs(fun, [file_out], kwargs)
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
            warn(f"{file_exc}, setting value to 'None'", category=QMFlows_Warning)
        else:
            self._results = results


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

    #: A class variable with the name of the generic .yaml file.
    #: Should be implemented by :class:`Package` subclasses.
    generic_dict_file: ClassVar[str] = NotImplemented

    #: A class variable with a special flag for used by the
    #: :class:`~qmflows.packages.package_wrapper.PackageWrapper` subclass.
    #: Used for denoting Job types without any generic .yaml files.
    generic_package: ClassVar[bool] = False

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

    @schedule(
        display="Running {self.pkg_name} {job_name}...",
        store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: Union[plams.Molecule, Chem.Mol],
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
        if self.generic_package:
            properties = package_properties[None]
        else:
            properties = package_properties[self.pkg_name]

        # Ensure that these variables have an actual value
        # Precaution against passing unbound variables to self.postrun()
        output_warnings = plams_mol = job_settings = None

        # There are not data from previous nodes in the dependecy trees
        # because of a failure upstream or the user provided None as argument
        if all(x is not None for x in [settings, mol]):
            #  Check if plams finishes normally
            try:
                # If molecule is an RDKIT molecule translate it to plams
                plams_mol = molkit.from_rdmol(mol) if isinstance(mol, Chem.Mol) else mol

                if job_name != '':
                    kwargs['job_name'] = job_name

                # Settings transformations
                self.prerun(settings, plams_mol, **kwargs)
                job_settings = self.generic2specific(settings, mol)

                # Run the job
                result = self.run_job(job_settings, plams_mol, **kwargs)

                # Check if there are warnings in the output that render the calculation
                # useless from the point of view of the user
                warnings_tolerance = kwargs.get("terminate_job_in_case_of_warnings")
                output_warnings = result.warnings

                if all(w is not None for w in [warnings_tolerance, output_warnings]):
                    issues = [w(msg) for msg, w in output_warnings.items()
                              if w in warnings_tolerance]
                    if issues:
                        warn(f"""
                        The Following Warning are rendered unacceptable in the Current
                        Workflow: {issues}\n
                        The results from Job: {job_name} are discarded.
                        """, category=QMFlows_Warning)
                        result = Result(None, None, job_name=job_name, dill_path=None,
                                        properties=properties, status='failed')

            # Otherwise pass an empty Result instance downstream
            except plams.core.errors.PlamsError as err:
                warn(f"Job {job_name} has failed.\n{err}", category=QMFlows_Warning)
                result = Result(None, None, job_name=job_name, dill_path=None,
                                properties=properties, status='failed')
        else:
            warn(f"""
            Job {job_name} has failed. Either the Settings or Molecule
            objects are None, probably due to a previous calculation failure
            """, category=QMFlows_Warning)

            # Send an empty object downstream
            result = Result(None, None, job_name=job_name, dill_path=None,
                            properties=properties, status='failed')

        # Label this calculation as failed if there are not dependecies coming
        # from upstream
        self.postrun(result, output_warnings, job_settings, plams_mol, **kwargs)
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
        generic_dict = self.get_generic_dict()

        specific_from_generic_settings = Settings()
        for k, v in settings.items():
            if k == "specific":
                continue
            elif k == 'input':  # Allow for PLAMS-style input; i.e. settings.input.blablabla
                specific_from_generic_settings.specific[self.pkg_name].update(v)
                continue

            key = generic_dict.get(k)
            if not key:
                self.handle_special_keywords(specific_from_generic_settings, k, v, mol)
                continue

            # The `key` variable can have four types of values:
            # * None
            # * A string
            # * A list with two elements; the second one being a dict
            # * A list of arbitrary length

            if isinstance(key, list):
                initial_key, *final_keys = key

                if isinstance(final_keys[0], dict):
                    v_tmp = final_keys[0].get(v)
                else:
                    v_tmp = Settings()
                    v_tmp.set_nested(final_keys, v)

                if v_tmp:
                    v = v_tmp
                key = initial_key

            if v:
                if isinstance(v, dict):
                    v = Settings(v)
                specific_from_generic_settings.specific[self.pkg_name][key] = v
            else:
                specific_from_generic_settings.specific[self.pkg_name][key] = Settings()

        return settings.overlay(specific_from_generic_settings)

    def get_generic_dict(self) -> Settings:
        """Load the .yaml file containing the translation from generic to the specific keywords of :attr:`Package.pkg_name`.

        Returns
        -------
        :class:`~qmflows.settings.Settings`
            A new Settings instance specific to :attr:`Package.pkg_name`.

        See Also
        --------
        :meth:`Package.generic2specific`
            Traverse *settings* and convert generic into package specific keys.

        """  # noqa
        try:
            path = join("data", "dictionaries", self.generic_dict_file)
        except TypeError as ex:
            if self.generic_dict_file is NotImplemented:
                raise NotImplementedError("The `Package.generic_dict_file` attribute should "
                                          "be implemented by `Package` subclasses") from ex
            raise ex

        str_yaml = pkg.resource_string("qmflows", path)
        return yaml2Settings(str_yaml)

    def __repr__(self) -> str:
        """Create a string representation of this instance.

        Returns
        -------
        :class:`str`
            A string representation of this instnce.

        """
        vars_str = ', '.join(f'{k}={v!r}' for k, v in sorted(vars(self).items()))
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

    @staticmethod
    @abstractmethod
    def run_job(settings: Settings, mol: plams.Molecule, job_name: str,
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
    \**kwargs : :data:`~typing.Any`
        Further keyword arguments to-be passed to :func:`call_default`.

    Returns
    -------
    :class:`Result`
        A new Result instance.

    See Also
    --------
    :func:`noodles.run_parallel`
        Run a workflow in parallel threads, storing results in a Sqlite3 database.

    """
    plams.init(path=path, folder=folder)
    plams.config.log.stdout = 0

    if runner is None:
        ret = call_default(job, kwargs.get('n_processes', 1),
                           kwargs.get('always_cache', True))
    else:
        raise ValueError(f"Don't know runner: {runner}")

    plams.finish()
    return ret


def call_default(wf: PromisedObject, n_processes: int, always_cache: bool) -> Result:
    """Run locally using several threads.

    Caching can be turned off by specifying ``cache=None``.

    """
    return run_parallel(
        wf, n_threads=n_processes, registry=registry,
        db_file='cache.db', always_cache=always_cache, echo_log=False)


class SerMolecule(Serialiser):
    """Based on the Plams molecule this class encode and decode the information related to the molecule using the JSON format."""  # noqa

    def __init__(self):
        super(SerMolecule, self).__init__(plams.Molecule)

    def encode(self, obj, make_rec):
        return make_rec(obj.as_dict())

    def decode(self, cls, data):
        return plams.Molecule.from_dict(data)


class SerMol(Serialiser):
    """Based on the RDKit molecule this class encodes and decodes the information related to the molecule using a string."""  # noqa

    def __init__(self):
        super(SerMol, self).__init__(Chem.Mol)

    def encode(self, obj, make_rec):
        return make_rec(base64.b64encode(obj.ToBinary()).decode('ascii'))

    def decode(self, cls, data):
        return Chem.Mol(base64.b64decode(data.encode('ascii')))


class SerSettings(Serialiser):
    """Class to encode and decode the :class:`~qmflows.Settings` class using its internal dictionary structure."""  # noqa

    def __init__(self):
        super(SerSettings, self).__init__(Settings)

    def encode(self, obj, make_rec):
        return make_rec(obj.as_dict())

    def decode(self, cls, data):
        return Settings(data)


#: A :class:`Registry` instance to-be returned by :func:`registry`.
REGISTRY: Registry = Registry(
    parent=serial.base() + arrays_to_hdf5(),
    types={
        Package: AsDict(Package),
        Path: SerPath(),
        plams.Molecule: SerMolecule(),
        Chem.Mol: SerMol(),
        Result: AsDict(Result),
        Settings: SerSettings(),
        plams.KFFile: SerReasonableObject(plams.KFFile),
        plams.KFReader: SerReasonableObject(plams.KFReader),
        np.floating: SerNumpyScalar(),
        np.integer: SerNumpyScalar()
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


def get_tmpfile_name() -> str:
    tmpfolder = join(plams.config.jm.workdir, 'tmpfiles')
    if not os.path.exists(tmpfolder):
        os.mkdir(tmpfolder)
    return join(tmpfolder, str(uuid.uuid4()))


def ignored_unused_kwargs(fun: Callable, args: Iterable, kwargs: Mapping) -> Any:
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
