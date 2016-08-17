# =======>  Standard and third party Python Libraries <======
from noodles import (Storable)
import pkg_resources as pkg
import plams

# ========================> Internal Modules  <================================
from qmworks.settings import Settings
from qmworks.packages.packages import Package, Result
from qmworks.fileFunctions import json2Settings
# ============================= Orca ==========================================


class _PropertyGetter(Storable):
    def __init__(self, result, section):
        super(_PropertyGetter, self).__init__()
        self.result = result
        self.section = section

    def __getattr__(self, prop):
        pass

    def __getitem__(self, prop):
        pass


class ORCA(Package):
    """
    This class prepare the input to run a Orca job using both Plams and
    templates. It also does the manangement of the input/output files resulting
    from running Orca and returns a Results object that containing the methods
    and data required to retrieve the output.
    """
    def __init__(self):
        super(ORCA, self).__init__("orca")
        self.generic_dict_file = 'generic2ORCA.json'

    def prerun(self):
        pass

    def run_job(self, settings, mol, job_name="ORCAjob"):

        orca_settings = Settings()
        orca_settings.input = settings.specific.orca
        print(orca_settings)
        result = plams.ORCAJob(molecule=mol, settings=orca_settings,
                               name=job_name).run()

        return ORCA_Result(orca_settings, mol, result, job_name)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        pass


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results"""

    def __init__(self, settings, molecule, result, job_name):
        self.settings = settings
        self._molecule = molecule
        self.result = result
        properties = 'data/dictionaries/propertiesORCA.json'
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        self.job_name = job_name

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.result.path
            "job_name": self.job_name}

    @classmethod
    def from_dict(cls, settings, molecule, filename, job_name):
        pass

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..

            dipole = result.properties.dipole

        """
        r = self.result.awk_output(script=self.prop_dict[prop])
        try:
            result = [float(i) for i in r]
        except:
            result = r
        if len(result) == 1:
            result = result[0]
        return result

    @property
    def molecule(self):
        pass

orca = ORCA()
