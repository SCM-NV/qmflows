
# =======>  Standard and third party Python Libraries <======
from warnings import warn
import plams

# ==================> Internal modules <====================
from qmworks.packages.packages import (Package, Result)
from qmworks.settings import Settings

# ==================> <======================
__all__ = ['dirac', 'DIRAC', 'DIRAC_Result']


class DIRAC(Package):
    """
    """
    def __init__(self):
        super().__init__("dirac")
        self.generic_dict_file = 'generic2Dirac.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, job_name="dirac_job"):

        dirac_settings = Settings()
        dirac_settings.input = settings.specific.dirac
        dirac_settings.ignore_molecule
        result = plams.DiracJob(name=job_name, settings=dirac_settings,
                                molecule=mol).run()

        return DIRAC_Result(dirac_settings, mol, result.job.name,
                            plams_dir=result.job.path)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Create the settings input for complex Dirac keywords
        """
        warn('Keyword ' + key + ' doesn\'t exist')

# Instance
dirac = DIRAC()


class DIRAC_Result(Result):
    """
    Class to access **DIRAC** Results.
    """
    def __init__(self, settings, molecule, job_name, plams_dir, project_name=None):
        properties = 'data/dictionaries/propertiesDIRAC.json'
        super().__init__(settings, molecule, job_name=job_name,
                         plams_dir=plams_dir, project_name=project_name,
                         properties=properties)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        plams_dir = archive["plams_dir"].path
        return DIRAC_Result(settings, molecule, job_name, plams_dir,
                            project_name)
