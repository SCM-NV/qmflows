
__all__ = ['gamess']

# =======>  Standard and third party Python Libraries <======
from qmworks.packages.packages import (Package, Result)
from qmworks.settings import Settings
from warnings import warn

import plams
# ======================================<>=====================================


class GAMESS(Package):
    """
    This class setup the requirement to run a Gamess-US Job
    <http://www.msg.ameslab.gov/gamess/>.
    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the ``rungms`` script.
    This class is not intended to be called directly by the user, instead the
    **gamess** function should be called.
    """
    def __init__(self):
        super(GAMESS, self).__init__("gamess")
        self.generic_dict_file = 'generic2gamess.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, work_dir=None, project_name=None,
                hdf5_file="quantum.hdf5", store_in_hdf5=True,
                job_name='gamess_job'):
        """
        Call the Cp2K binary using plams interface.

        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the
        numerical results.
        :type hdf5_file: String
        :param input_file_name: Optional name for the input.
        :type input_file_name: String
        :param out_file_name: Optional name for the output.
        :type out_file_name: String
        :param store_in_hdf5: wether to store the output arrays in HDF5 format.
        :type store_in_hdf5: Bool
        """
        gamess_settings = Settings()
        gamess_settings.input = settings.specific.gamess
        job = plams.GamessJob(molecule=mol, name=job_name,
                              settings=gamess_settings)
        r = job.run()

        result = Gamess_Result(gamess_settings, mol, r.job.name,
                               plams_dir=r.job.path, work_dir=work_dir,
                               path_hdf5=hdf5_file, project_name=project_name,
                               status=job.status)

        return result

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Create the settings input for complex cp2k keys

        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param key: Special key declared in ``settings``.
        :param value: Value store in ``settings``.
        :param mol: molecular Geometry
        :type mol: plams Molecule
        """
        warn('Keyword ' + key + ' doesn\'t exist')


class Gamess_Result(Result):
    """
    Class providing access to CP2K result.
    """
    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 work_dir=None, path_hdf5=None, project_name=None,
                 properties='data/dictionaries/propertiesGAMESS.json',
                 status='done'):
        super().__init__(settings, molecule, job_name, plams_dir,
                         work_dir=work_dir, path_hdf5=path_hdf5,
                         project_name=None, properties=properties,
                         status=status)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        """
        Create a :class:`~CP2K_Result` instance using the data serialized in
        a dictionary.

        :param cls:
        :param settings: Job Settings.
        :param molecule: molecular Geometry.
        :param job_name: Name of the job.
        :param plams_dir: Absolute path to plams output folder
        :param archive: dictionary containing the paths to the input/output
        folders.
        :param path_hdf5: Path to the HDF5 file that contains the numerical
        results.
        """
        plams_dir = archive.get("plams_dir").path
        work_dir = archive.get("work_dir")
        return Gamess_Result(settings, molecule, job_name,
                             plams_dir=plams_dir, work_dir=work_dir,
                             project_name=project_name)

    def __getattr__(self, prop):
        """Returns a section of the results.
        The property is extracted from a  result file, which is recursively
        search for in the GAMESS settings

        Example:
        ..
            Hessian_matrix = result.hessian
        """
        result = super().__getattr__(prop)
        if result is None:
            msg = """
            Maybe you need to provided to the gamess
            function the optional keyword 'work_dir' containing the path
            to the SCR folder where GAMESS stores the *.dat and other
            output files"""
            warn(msg)

        return result


gamess = GAMESS()
