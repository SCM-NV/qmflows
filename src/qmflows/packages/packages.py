"""Common funcionality to call all the quantum packages."""

import base64
import fnmatch
import importlib
import inspect
import os
import uuid
import warnings
from types import ModuleType
from pathlib import Path
from functools import partial
from os.path import join
from warnings import warn
from typing import (Any, Callable, Optional, Dict, Union,
                    Iterable, Mapping, Iterator, Type)

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

from ..fileFunctions import json2Settings
from ..settings import Settings

__all__ = ['package_properties',
           'Package', 'run', 'registry', 'Result',
           'SerMolecule', 'SerSettings']

base_path = Path('data') / 'dictionaries'

#: A dictionary mapping package names to .json files.
package_properties: Dict[str, Path] = {
    'adf': base_path / 'propertiesADF.json',
    'dftb': base_path / 'propertiesDFTB.json',
    'cp2k': base_path / 'propertiesCP2K.json',
    'dirac': base_path / 'propertiesDIRAC.json',
    'gamess': base_path / 'propertiesGAMESS.json',
    'orca': base_path / 'propertiesORCA.json'
}
del base_path

WarnMap = Mapping[str, Type[Warning]]
WarnDict = Dict[str, Type[Warning]]
ParserFunc = Callable[[str, WarnMap], Optional[WarnDict]]


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
        :param properties: path to the `JSON` file containing the data to
                           load the parser on the fly.
        :type properties: str

        """
        plams_dir = None if plams_dir is None else Path(plams_dir)
        self.settings = settings
        self._molecule = molecule
        xs = pkg.resource_string("qmflows", str(properties))
        self.prop_dict = json2Settings(xs)
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
                warn(f"Generic property {prop!r} not defined")

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
            """)
        return None

    def get_property(self, prop: str) -> Any:
        """Look for the optional arguments to parse a property, which are stored in the properties dictionary."""  # noqa
        # Read the JSON dictionary than contains the parsers names
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
    def results(self) -> plams.Results:
        """Getter for :attr:`Result.results`.

        Get will load the .dill file and add all of its class attributes to this instance,
        barring the following three exceptions:

        * Private attributes/methods.
        * Magic methods.
        * Methods/attributes with names already contained within this instance.

        """
        if self._results_open:
            return self._results
        else:
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
            warnings.simplefilter("ignore", category=UserWarning)

            # Unpickle the results
            file_err = None
            try:
                results = plams.load(self._results).results
                assert results is not None, f'Failed to unpickle {self._results}'
            except (AssertionError, plams.FileError) as ex:
                file_err = ex
            else:
                attr_set = set(dir(self))
                for name in dir(results):
                    if name.startswith('_') or name in attr_set:
                        continue  # Skip methods which are either private, magic or preexisting

                    results_func = getattr(results, name)
                    setattr(self, name, results_func)

        # Failed to find or unpickle the .dill file; issue a warning
        if file_err is not None:
            self._results = None
            warn(str(file_err))
        else:
            self._results = results


@has_scheduled_methods
class Package:
    """|Package| is the base class to handle the invocation to different quantum package.

    The only relevant attribute of this class is :attr:`Package.pkg_name` which is a
    string representing the quantum package name that is going to be used to
    carry out the compuation.

    """

    def __init__(self, pkg_name: str) -> None:
        """Initialize a :class:`Package` instance."""
        self.pkg_name = pkg_name
        self.generic_dict_file = None  # will raise a NotImplementedError if .generic_dict_file

    @schedule(
        display="Running {self.pkg_name} {job_name}...",
        store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: Union[plams.Molecule, Chem.Mol],
                 job_name: str = '', **kwargs: Any) -> Result:
        """Perform a job with the package specified by :attr:`Package.pkg_name`.

        :parameter settings: user settings
        :type settings: |Settings|
        :parameter mol: Molecule to run the calculation.
        :type mol: plams Molecule

        """
        properties = package_properties[self.pkg_name]

        # There are not data from previous nodes in the dependecy trees
        # because of a failure upstream or the user provided None as argument
        if all(x is not None for x in [settings, mol]):
            #  Check if plams finishes normally
            try:
                # If molecule is an RDKIT molecule translate it to plams
                mol = molkit.from_rdmol(mol) if isinstance(mol, Chem.Mol) else mol

                if job_name != '':
                    kwargs['job_name'] = job_name

                # Settings transformations
                job_settings = self.generic2specific(settings, mol)

                # Run the job
                self.prerun()
                result = self.run_job(job_settings, mol, **kwargs)

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
                        """)
                        result = Result(None, None, job_name=job_name, dill_path=None,
                                        properties=properties, status='failed')

            # Otherwise pass an empty Result instance downstream
            except plams.core.errors.PlamsError as err:
                warn(f"Job {job_name} has failed.\n{err}")
                result = Result(None, None, job_name=job_name, dill_path=None,
                                properties=properties, status='failed')
        else:
            warn(f"""
            Job {job_name} has failed. Either the Settings or Molecule
            objects are None, probably due to a previous calculation failure
            """)

            # Send an empty object downstream
            result = Result(None, None, job_name=job_name, dill_path=None,
                            properties=properties, status='failed')

        # Label this calculation as failed if there are not dependecies coming
        # from upstream
        self.postrun()
        return result

    def generic2specific(self, settings: Settings,
                         mol: Optional[plams.Molecule] = None) -> Settings:
        """Traverse ``settings`` and convert generic into package specific keys.

        Traverse all the key, value pairs of the ``settings``, translating
        the generic keys into package specific keys as defined in the specific
        dictionary. If one key is not in the specific dictionary an error
        is raised. These new specific settings take preference over existing
        specific settings.

        :parameter settings: Settings provided by the user.
        :type      settings: Settings
        :parameter mol: Molecule to run the calculation.
        :type mol: plams Molecule

        """
        generic_dict = self.get_generic_dict()

        specific_from_generic_settings = Settings()
        for k, v in settings.items():
            if k != "specific":
                key = generic_dict.get(k)
                if key:
                    if isinstance(key, list):
                        if isinstance(key[1], dict):
                            value = key[1][v]
                        else:
                            value = {key[1]: v}
                        if value:
                            v = value
                        key = key[0]
                    if v:
                        if isinstance(v, dict):
                            v = Settings(v)
                        specific_from_generic_settings \
                            .specific[self.pkg_name][key] = v
                    else:
                        specific_from_generic_settings \
                            .specific[self.pkg_name][key]
                else:
                    self.handle_special_keywords(
                        specific_from_generic_settings, k, v, mol)
        return settings.overlay(specific_from_generic_settings)

    def get_generic_dict(self) -> Settings:
        """Load the JSON file containing the translation from generic to the specific keywords of :attr:`Pacakge.self.pkg_name``."""  # noqa
        try:
            path = join("data", "dictionaries", self.generic_dict_file)
        except TypeError as ex:
            if self.generic_dict_file is None:
                raise NotImplementedError("The `Package.generic_dict_file` attribute "
                                          "should be implemented by Package subclasses")
            raise ex

        str_json = pkg.resource_string("qmflows", path)
        return json2Settings(str_json)

    def __repr__(self) -> str:
        """Return a :class:`str` representation of this instance."""
        vars_str = ', '.join(f'{k}={v!r}' for k, v in sorted(vars(self).items()))
        return f'{self.__class__.__name__}({vars_str})'

    def prerun(self) -> None:
        """Run a set of tasks before running the actual job."""
        pass

    def postrun(self) -> None:
        """Run a set of tasks after running the actual job."""
        pass

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Abstract method; should be implemented by the child class."""
        raise NotImplementedError("trying to call an abstract method")

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_name: Optional[str] = None,
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> Result:
        """Abstract method; should be implemented by the child class."""
        raise NotImplementedError("The class representing a given quantum packages "
                                  "should implement this method")


def run(job: PromisedObject, runner: Optional[str] = None,
        path: Union[None, str, os.PathLike] = None,
        folder: Union[None, str, os.PathLike] = None,
        **kwargs: Any) -> Result:
    """Pickup a runner and initialize it.

    :params job: computation to run
    :type job: Promise Object
    :param runner: Type of runner to use
    :type runner: String

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
    defaults = list(filter(lambda t: t[1].default != inspect._empty,
                           ps.items()))
    # there are not keyword arguments in the function
    if not kwargs or not defaults:
        return fun(*args)
    else:  # extract from kwargs only the used keyword arguments
        d = {k: kwargs[k] for k, _ in defaults}
        return fun(*args, **d)


def parse_output_warnings(job_name: str,
                          plams_dir: Union[None, str, os.PathLike],
                          parser: ParserFunc,
                          package_warnings: WarnMap) -> Optional[WarnDict]:
    """Look out for warnings in the output file."""
    output_files = find_file_pattern('*out', plams_dir)
    try:
        return parser(next(output_files), package_warnings)
    except StopIteration:
        warn(f"job: {job_name} has failed. check folder: {plams_dir}")
        return None
