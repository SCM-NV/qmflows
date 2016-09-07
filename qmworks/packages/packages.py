
# ========>  Standard and third party Python Libraries <======
from os.path import join
from rdkit import Chem
import subprocess
import base64

import os
import plams
import pkg_resources as pkg

# ==================> Internal modules <====================
from noodles import (schedule_hint, has_scheduled_methods, serial)
from noodles.display import (NCDisplay)
from noodles.files.path import (Path, SerPath)
from noodles.run.run_with_prov import run_parallel_opt
from noodles.serial import (Serialiser, Registry, AsDict)
from noodles.serial.base import SerAutoStorable

from qmworks.settings import Settings
from qmworks import rdkitTools
from qmworks.fileFunctions import json2Settings

# ==============================================================
__all__ = ['Package', 'run', 'registry', 'Result',
           'SerMolecule', 'SerSettings']


class Result:

    def __init__(self, settings, molecule, job_name, plams_dir, work_dir=None,
                 path_hdf5=None, project_name=None, properties=None):
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
        :param hdf5_file: path to the file containing the numerical results.
        :type hdf5_file: str
        :param properties: path to the `JSON` file containing the properties
        addresses inside the `HDF5` file.
        :type properties: str
        """
        self.settings = settings
        self._molecule = molecule
        self.hdf5_file = path_hdf5
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        self.archive = {"plams_dir": Path(plams_dir),
                        'work_dir': work_dir}
        self.project_name = project_name
        self.job_name = job_name

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.archive,
            "job_name": self.job_name}

    def awk_output(self, script='', progfile=None, **kwargs):
        """awk_output(script='', progfile=None, **kwargs)
        Shortcut for :meth:`~Results.awk_file` on the output file."""
        output = self.name + ".out"
        return self.awk_file(output, script, progfile, **kwargs)

    def awk_file(self, filename, script='', progfile=None, **kwargs):
        """awk_file(filename, script='', progfile=None, **kwargs)
        Execute an AWK script on a file given by *filename*.

        The AWK script can be supplied in two ways: either by directly passing
        the contents of the script (should be a single string) as a *script*
        argument, or by providing the path (absolute or relative to the file
        pointed by *filename*) to some external file containing the actual AWK
        script using *progfile* argument. If *progfile* is not ``None``, the
        *script* argument is ignored.

        Other keyword arguments (*\*\*kwargs*) can be used to pass additional
        variables to AWK (see ``-v`` flag in AWK manual)

        Returned value is a list of lines (strings). See ``man awk`` for details.
        """
        cmd = ['awk']
        for k, v in kwargs.items():
            cmd += ['-v', '%s=%s' % (k, v)]
        if progfile:
            if os.path.isfile(progfile):
                cmd += ['-f', progfile]
            else:
                raise FileError('File %s not present' % progfile)
        else:
            cmd += [script]
        ret = subprocess.check_output(cmd + [filename], cwd=self.path).decode('utf-8').split('\n')
        if ret[-1] == '':
            ret = ret[:-1]
        try:
            result = [float(i) for i in ret]
        except:
            result = ret
        if len(result) == 1:
            result = result[0]
        return result


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

        if isinstance(mol, Chem.Mol):
            mol = rdkitTools.rdkit2plams(mol)

        if job_name != '':
            kwargs['job_name'] = job_name

        self.prerun()

        job_settings = self.generic2specific(settings, mol)

        result = self.run_job(job_settings, mol, **kwargs)

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


def run(job, runner=None, **kwargs):
    """
    Pickup a runner and initialize it.

    :params job: computation to run
    :type job: Promise Object
    :param runner: Type of runner to use
    :type runner: String
    """

    if runner is None:
        return call_default(job, **kwargs)
    elif runner.lower() is 'xenon':
        return call_xenon(job, **kwargs)


def call_default(job, n_processes=1):
    """
    Run locally using several threads.
    """
    with NCDisplay() as display:
        return run_parallel_opt(
            job, n_threads=n_processes,
            registry=registry, jobdb_file='cache.json',
            display=display)


def call_xenon(job, **kwargs):
    """
    See :
        https://github.com/NLeSC/Xenon-examples/raw/master/doc/tutorial/xenon-tutorial.pdf
    """
    pass
    # nproc = kwargs.get('n_processes')
    # nproc = nproc if nproc is not None else 1

    # xenon_config = XenonConfig(jobs_scheme='local')

    # job_config = RemoteJobConfig(registry=serial.base, time_out=1)

    # return run_xenon(job, nproc, xenon_config, job_config)


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
        parent=serial.base(),
        types={
            Package: AsDict(Package),
            Path: SerPath(),
            plams.Molecule: SerMolecule(),
            Chem.Mol: SerMol(),
            Result: SerAutoStorable(Result),
            Settings: SerSettings()})
