__all__ = ['dirac']

from qmflows.packages.packages import (Package, package_properties, Result)
from qmflows.settings import Settings
from scm import plams
from warnings import warn


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
        job = plams.interfaces.dirac.DiracJob(name=job_name,
                                              settings=dirac_settings,
                                              molecule=mol)
        result = job.run()

        # Relative job path
        relative_plams_path = '/'.join(result.job.path.split('/')[-2:])

        return DIRAC_Result(dirac_settings, mol, result.job.name,
                            plams_dir=relative_plams_path, status=job.status)

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
    def __init__(self, settings, molecule, job_name, plams_dir,
                 status='done', warnings=None):
        properties = package_properties['dirac']
        super().__init__(settings, molecule, job_name=job_name,
                         plams_dir=plams_dir, properties=properties,
                         status=status, warnings=warnings)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status, warnings):
        plams_dir = archive["plams_dir"].path
        return DIRAC_Result(settings, molecule, job_name, plams_dir,
                            status, warnings)
