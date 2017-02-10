
# ========>  Standard and third party Python Libraries <======
from functools import partial
from os.path import join
from rdkit import Chem
from typing import (Any, Callable, Dict, List)

import base64
import fnmatch
import importlib
import inspect
import os
import plams
import pkg_resources as pkg
import builtins

# ==================> Internal modules <====================
from noodles import (schedule_hint, has_scheduled_methods, serial)
from noodles.display import (NCDisplay)
from noodles.files.path import (Path, SerPath)
from noodles.run.run_with_prov import run_parallel_opt
from noodles.serial import (Serialiser, Registry, AsDict)
from noodles.serial.base import SerStorable
from noodles.run.xenon import (
    XenonKeeper, XenonConfig, RemoteJobConfig, run_xenon_prov)
from noodles.serial.numpy import arrays_to_hdf5

from qmworks.settings import Settings
from qmworks import molkit
from qmworks.fileFunctions import json2Settings
from qmworks.utils import concatMap
from warnings import warn
# ==============================================================
__all__ = ['import_parser', 'package_properties',
           'Package', 'run', 'registry', 'Result',
           'SerMolecule', 'SerSettings']

package_properties = {
    'adf': 'data/dictionaries/propertiesADF.json',
    'dftb': 'data/dictionaries/propertiesDFTB.json',
    'cp2k': 'data/dictionaries/propertiesCP2K.json',
    'dirac': 'data/dictionaries/propertiesDIRAC.json',
    'gamess': 'data/dictionaries/propertiesGAMESS.json',
    'orca': 'data/dictionaries/propertiesORCA.json'
}


class Result:
    """
    Class containing the result associated with a quantum chemistry simulation.
    """
    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 work_dir=None, properties=None, status='done'):
        """
        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param molecule: molecular Geometry
        :type molecule: plams Molecule
        :param job_name: Name of the computations
        :type job_name: str
        :param plams_dir: path to the ``Plams`` folder.
        :type plams_dir: str
        :param work_dir: scratch or another directory different from
        the `plams_dir`.
        type work_dir: str
        :param properties: path to the `JSON` file containing the data to
                           load the parser on the fly.
        :type properties: str
        """
        self.settings = settings
        self._molecule = molecule
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        self.archive = {"plams_dir": Path(plams_dir),
                        'work_dir': work_dir}
        self.job_name = job_name
        self.status = status

    def as_dict(self):
        """
        Method to serialize as a JSON dictionary the results given
        by an ``Package`` computation.
        """
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "job_name": self.job_name,
            "archive": self.archive,
            "status": self.status}

    def from_dict(cls, settings, molecule, job_name, archive,
                  status):
        """
        Methods to deserialize an `Result`` object.
        """
        raise NotImplementedError()

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..
            dipole = result.dipole
        """
        crash_status = ['failed', 'crashed']
        is_private = prop.startswith('__') and prop.endswith('__')
        # if self.status == 'successful':
        if self.status not in crash_status and prop in self.prop_dict:
            return self.get_property(prop)
        elif (self.status not in crash_status and not is_private and
              prop not in self.prop_dict):
            msg = "Generic property '" + str(prop) + "' not defined"
            warn(msg)
            return None
        elif (self.status in crash_status) and not is_private:
            warn("""
            It is not possible to retrieve property: '{}'
            Because Job: '{}' has failed. Check the output.\n
            Are you sure that you have the package installed or
             you have loaded the package in the cluster. For example:
            `module load AwesomeQuantumPackage/3.1421`
            """.format(prop, self.job_name))
            return None

    def get_property(self, prop):
        """
        Look for the optional arguments to parse a property, which are stored
        in the properties dictionary.
        """
        # Read the JSON dictionary than contains the parsers names
        ds = self.prop_dict[prop]

        # extension of the output file containing the property value
        file_ext = ds['file_ext']

        # If there is not work_dir returns None
        work_dir = self.archive.get('work_dir')

        # Plams dir
        plams_dir = self.archive['plams_dir'].path

        # Search for the specified output file in the folders
        file_pattern = ds.get('file_pattern')
        if file_pattern is None:
            file_pattern = '{}*.{}'.format(self.job_name, file_ext)

        output_files = concatMap(partial(find_file_pattern, file_pattern),
                                 [plams_dir, work_dir])
        if output_files:
            file_out = output_files[0]
            fun = getattr(import_parser(ds), ds['function'])
            # Read the keywords arguments from the properties dictionary
            kwargs = ds.get('kwargs') if ds.get('kwargs') is not None else {}
            kwargs['plams_dir'] = plams_dir
            return ignored_unused_kwargs(fun, [file_out], kwargs)
        else:
            msg = """
            Property {} not found. No output file called: {}. Folder used:
            plams_dir = {}\n
            work_dir {}\n
            """.format(prop, file_pattern, plams_dir, work_dir)
            raise FileNotFoundError(msg)


@has_scheduled_methods
class Package:
    """
    |Package| is the base class to handle the invocation to different
    quantum package.
    The only relevant attribute of this class is ``self.pkg_name`` which is a
    string representing the quantum package name that is going to be used to
    carry out the compuation.

    Only two arguments are required
    """
    def __init__(self, pkg_name):
        super(Package, self).__init__()
        self.pkg_name = pkg_name

    @schedule_hint(
        display="Running {self.pkg_name} {job_name}...",
        store=True, confirm=True)
    def __call__(self, settings, mol, job_name='', **kwargs):
        """
        This function performs a job with the package specified by
        self.pkg_name

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
                if isinstance(mol, Chem.Mol):
                    mol = molkit.from_rdmol(mol)

                if job_name != '':
                    kwargs['job_name'] = job_name

                # Settings transformations
                job_settings = self.generic2specific(settings, mol)

                # Run the job
                self.prerun()
                result = self.run_job(job_settings, mol, **kwargs)
            # Otherwise pass an empty Result instance
            except plams.PlamsError:
                result = Result(None, None, job_name=job_name,
                                properties=properties, status='failed')
            # except Exception as e:
            #     print("Exception e: ", type(e), e.args)
        else:
            result = Result(None, None, job_name=job_name,
                            properties=properties, status='failed')

        # Label this calculation as failed if there are not dependecies coming
        # from upstream

        self.postrun()

        return result

    def generic2specific(self, settings, mol=None):
        """
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
                            value = key[1]
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

    def get_generic_dict(self):
        """
        Loads the JSON file containing the translation from generic to
        the specific keywords of ``self.pkg_name``.
        """
        path = join("data/dictionaries", self.generic_dict_file)
        str_json = pkg.resource_string("qmworks", path)

        return json2Settings(str_json)

    def __str__(self):
        return self.pkg_name

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        This method should be implemented by the child class.
        """
        msg = "trying to call an abstract method"
        raise NotImplementedError(msg)

    @staticmethod
    def run_job(settings, mol, job_name=None, work_dir=None, **kwargs):
        """
        This method should be implemented by the child class.
        """
        msg = "The class representing a given quantum packages should \
        implement this method"
        raise NotImplementedError(msg)


def run(job, runner=None, path=None, folder=None, **kwargs):
    """
    Pickup a runner and initialize it.

    :params job: computation to run
    :type job: Promise Object
    :param runner: Type of runner to use
    :type runner: String
    """

    initialize = False
    try:
        builtins.config
        if path and os.path.abspath(path) != builtins.config.jm.path or \
                folder and folder != builtins.config.jm.folder:
            msg = "Reinitializing Plams with new path and/or folder name.\n"
            warn(msg)
            plams.finish()
            plams.init(path=path, folder=folder)
    except:
        plams.init(path=path, folder=folder)
        initialize = True
    builtins.config.log.stdout = 0
    builtins.config.jobmanager.jobfolder_exists = 'rename'
    if runner is None:
        ret = call_default(job, **kwargs)
    elif runner.lower() == 'xenon':
        ret = call_xenon(job, **kwargs)
    else:
        raise "Don't know runner: {}".format(runner)
    if initialize:
        plams.finish()
    return ret


def call_default(job, n_processes=1, cache='cache.json'):
    """
    Run locally using several threads.
    """
    with NCDisplay() as display:
        return run_parallel_opt(
            job, n_threads=n_processes,
            registry=registry, jobdb_file=cache,
            display=display)


def call_xenon(job, n_processes=1, cache='cache.json', user_name=None, adapter='slurm',
               queue_name=None, host_name=None, workdir=None, timeout=60000, **kwargs):
    """
    See :
        https://github.com/NLeSC/Xenon-examples/raw/master/doc/tutorial/xenon-tutorial.pdf
    """
    dict_properties = {
        'slurm': {'xenon.adaptors.slurm.ignore.version': 'true'},
        'pbs': {'xenon.adaptors.pbs.ignore.version': 'true'}
    }
    with XenonKeeper(log_level='DEBUG') as Xe:
        certificate = Xe.credentials.newCertificateCredential(
            'ssh', os.environ["HOME"] + '/.ssh/id_rsa', user_name, '', None)

        xenon_config = XenonConfig(
            jobs_scheme=adapter,
            location=host_name,
            credential=certificate,
            jobs_properties=dict_properties[adapter]
        )
        print(xenon_config.__dict__)

        if workdir is None:
            workdir = '/home/' + user_name

        job_config = RemoteJobConfig(
            registry=registry,
            init=plams.init,
            finish=plams.finish,
            queue=queue_name,
            time_out=timeout,
            working_dir=workdir
        )

        with NCDisplay() as display:
            result = run_xenon_prov(
                job, Xe, cache, n_processes,
                xenon_config, job_config, display=display)

    return result


class SerMolecule(Serialiser):
    """
    Based on the Plams molecule this class encode and decode the
    information related to the molecule using the JSON format.
    """
    def __init__(self):
        super(SerMolecule, self).__init__(plams.Molecule)

    def encode(self, obj, make_rec):
        return make_rec(obj.as_dict())

    def decode(self, cls, data):
        return plams.Molecule.from_dict(**data)


class SerMol(Serialiser):
    """
    Based on the RDKit molecule this class encodes and decodes the
    information related to the molecule using a string.
    """
    def __init__(self):
        super(SerMol, self).__init__(Chem.Mol)

    def encode(self, obj, make_rec):
        return make_rec(base64.b64encode(obj.ToBinary()).decode('ascii'))

    def decode(self, cls, data):
        return Chem.Mol(base64.b64decode(data.encode('ascii')))


class SerSettings(Serialiser):
    """
    Class to encode and decode the ~qmworks.Settings class using
    its internal dictionary structure.
    """

    def __init__(self):
        super(SerSettings, self).__init__(Settings)

    def encode(self, obj, make_rec):
        return make_rec(obj.as_dict())

    def decode(self, cls, data):
        return Settings(data)


def registry():
    """
    This function pass to the noodles infrascture all the information
    related to the Structure of the Package object that is schedule.
    This *Registry* class contains hints that help Noodles to encode
    and decode this Package object.
    """
    return Registry(
        parent=serial.base() + arrays_to_hdf5(),
        types={
            Package: AsDict(Package),
            Path: SerPath(),
            plams.Molecule: SerMolecule(),
            Chem.Mol: SerMol(),
            Result: SerStorable(Result),
            Settings: SerSettings()})


def import_parser(ds, module_root="qmworks.parsers"):
    """
    Import parser for the corresponding property.
    """
    module_sufix = ds['parser']
    module_name = module_root + '.' + module_sufix

    return importlib.import_module(module_name)


def find_file_pattern(pat, folder):
    if folder is not None and os.path.exists(folder):
        return map(lambda x: join(folder, x),
                   fnmatch.filter(os.listdir(folder), pat))
    else:
        return []


def ignored_unused_kwargs(fun: Callable, args: List, kwargs: Dict) -> Any:
    """
    Inspect the signature of function `fun` and filter the keyword arguments,
    which are the ones that have a nonempty default value. Then extract
    from the dict `kwargs` those key-value pairs ignoring the rest.
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
