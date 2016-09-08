# =======>  Standard and third party Python Libraries <======
from noodles import (Storable)
from qmworks.settings import Settings
from qmworks.packages.packages import Package, Result

import plams
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

        result = plams.ORCAJob(molecule=mol, settings=orca_settings,
                               name=job_name).run()

        return ORCA_Result(orca_settings, mol, result.job.path, result.job.name)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        pass


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results"""

    def __init__(self, settings, molecule, job_name, plams_dir, work_dir=None,
                 path_hdf5=None, project_name=None,
                 properties='data/dictionaries/propertiesORCA.json'):
        super().__init__(settings, molecule, job_name, plams_dir,
                         work_dir=work_dir, path_hdf5=path_hdf5,
                         project_name=project_name, properties=properties)

    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "archive": self.path,
            "name": self.name}
    
    @classmethod
    def from_dict(cls, settings, molecule, archive, job_name):
        return ORCA_Result(settings, molecule, archive, job_name)

    def __getattr__(self, prop):
        """Returns a section of the results.

        Example:

        ..

            dipole = result.dipole

        """
        return self.awk_output(script=self.prop_dict[prop])

    @property
    def molecule(self):
        pass

orca = ORCA()
