
__all__ = ['Package', 'run', 'registry', 'Result',
           'SerMolecule', 'SerSettings']

# ========>  Standard and third party Python Libraries <======
from os.path import join

import plams
import pkg_resources as pkg

# ==================> Internal modules <====================
from noodles import (schedule_hint, has_scheduled_methods, serial,
                     run_parallel)
from noodles.serial import (Serialiser, Registry, AsDict)
from noodles.serial.base import SerAutoStorable

from qmworks.settings import Settings
from qmworks.fileFunctions import json2Settings

# ==============================================================


class Result:
    def __init__(self):
        pass


@has_scheduled_methods
class Package:
    """
    |Package| is the base class to handle the invocation to different
    quantum package.
    The only relevant attribute of this class is ``self.pkg_name`` which is a
    string representing the quantum package name that is going to be used to carry
    out the compuation.
    Only two arguments are required
    """
    def __init__(self, pkg_name):
        super(Package, self).__init__()
        self.pkg_name = pkg_name

    @schedule_hint(display="Running {self.pkg_name} ...")
    def __call__(self, settings, mol, **kwargs):
        """
        This function performs a job with the package specified by self.pkg_name

        :parameter settings: user settings
        :type settings: |Settings|
        :parameter mol: Molecule to run the calculation.
        :type mol: plams Molecule
        """

        self.prerun()

        job_settings = self.generic2specific(settings, mol)
        result = self.run_job(job_settings, mol, **kwargs)

        self.postrun()

        return result

    def generic2specific(self, settings, mol=None):
        """
        Traverse all the key, value pairs of the ``settings``, translating
        the generic keys into package specific keys as defined in the specific
        dictionary. If one key is not in the specific dictionary an error is raised
        These new specific settings take preference over existing specific settings.

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
                        specific_from_generic_settings.specific[self.pkg_name][key] = v
                    else:
                        specific_from_generic_settings.specific[self.pkg_name][key]
                else:
                    self.handle_special_keywords(specific_from_generic_settings, k, v, mol)
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
    return run_parallel(job, n_threads=n_processes)


def call_xenon(job, **kwargs):
    """
    See : https://github.com/NLeSC/Xenon-examples/raw/master/doc/tutorial/xenon-tutorial.pdf
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
            plams.Molecule: SerMolecule(),
            Result: SerAutoStorable(Result),
            Settings: SerSettings()})

