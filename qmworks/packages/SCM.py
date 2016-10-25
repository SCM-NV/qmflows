
# =======>  Standard and third party Python Libraries <======

from os.path import join
from warnings import warn
from qmworks.settings import Settings
from qmworks.packages.packages import Package, Result  # ChemResult

import builtins
import plams
# ========================= ADF ============================


class ADF(Package):
    """
    This class takes care of calling the *ADF* quantum package.
    it uses both Plams and the Templates module to create the input
    files, while Plams takes care of running the Job.
    It returns a ``ADF_Result`` object containing the output data.
    """
    def __init__(self):
        super(ADF, self).__init__("adf")
        self.generic_dict_file = 'generic2ADF.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, job_name='ADFjob'):
        """
        Execute ADF job.

        :param settings: user input settings.
        :type settings: |Settings|
        :param mol: Molecule to run the simulation
        :type mol: Plams Molecule
        :parameter input_file_name: The user can provide a name for the
                                    job input.
        :type input_file_name: String
        :parameter out_file_name: The user can provide a name for the
                                  job output.
        :type out_file_name: String
        :returns: :class:`~qmworks.packages.SCM.ADF_Result`
        """
        adf_settings = Settings()
        adf_settings.input = settings.specific.adf
        result = plams.ADFJob(name=job_name, molecule=mol,
                              settings=adf_settings).run()
        path_t21 = result._kf.path

        return ADF_Result(adf_settings, mol, result.job.name, path_t21,
                          plams_dir=result.job.path)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        some keywords provided by the user do not have a straightforward
        translation to *ADF* input and require some hooks that handles the
        special behaviour of the following keywords:

        * ``freeze``
        * ``selected_atoms``
        """
        if key == "freeze":
            settings.specific.adf.geometry.optim = "cartesian"
            for a in value:
                settings.specific.adf.constraints['atom ' + str(a + 1)] = ""
        elif key == "selected_atoms":
            settings.specific.adf.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in range(len(mol)):
                    if a not in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol.atoms[a].symbol not in value:
                        name = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[name] = ""
        elif key == "inithess":
            hess_path = builtins.config.jm.workdir + "/tmp_hessian.txt"
            hess_file = open(hess_path, "w")
            hess_file.write(" ".join(['{:.6f}'.format(v) for v in value]))
            settings.specific.adf.geometry.inithess = hess_path
        else:
            msg = 'Generic keyword "' + key + '" not implemented for package ADF.'
            warn(msg)


class ADF_Result(Result):
    """Class providing access to PLAMS ADFJob result results"""

    def __init__(self, settings, molecule, job_name, path_t21, plams_dir=None,
                 project_name=None):
        properties = 'data/dictionaries/propertiesADF.json'
        super().__init__(settings, molecule, job_name,
                         plams_dir=plams_dir, project_name=project_name,
                         properties=properties)
        self.kf = plams.kftools.KFFile(path_t21)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        """
        Methods to deserialize an `ADF_Result` object.
        """
        plams_dir = archive["plams_dir"].path
        path_t21 = join(plams_dir, '{}.t21'.format(job_name))
        return ADF_Result(settings, molecule, job_name, path_t21, plams_dir,
                          project_name)

    def get_property_kf(self, prop, section=None):
        return self.kf.read(section, prop)

    @property
    def molecule(self, unit='bohr', internal=False, n=1):
        """WARNING: Cheap copy from PLAMS, do not keep this!!!"""
        m = self._molecule.copy()
        natoms = len(m)
        # Find out correct location
        coords = self.kf.read('Geometry', 'xyz InputOrder')
        coords = [coords[i: i + 3] for i in range(0, len(coords), 3)]

        if len(coords) > natoms:
            coords = coords[(n - 1) * natoms: n * natoms]
        if internal:
            mapping = self._int2inp()
            coords = [coords[mapping[i] - 1] for i in range(len(coords))]
        for at, coord in zip(m, coords):
            at.move_to(coord, unit)

        return m


class DFTB(Package):
    """
    Add some documentation to this class
    """
    def __init__(self):
        super().__init__("dftb")
        self.generic_dict_file = 'generic2DFTB.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, job_name='DFTBjob'):
        """
        Execute an DFTB job with the *ADF* quantum package.

        :param settings: user input settings.
        :type settings: |Settings|
        :param mol: Molecule to run the simulation
        :type mol: Plams Molecule
        :parameter input_file_name: The user can provide a name for the
                                   job input.
        :type input_file_name: String
        :parameter out_file_name: The user can provide a name for the
                                 job output.
        :type out_file_name: String
        :returns: :class:`~qmworks.packages.SCM.DFTB_Result`
        """
        dftb_settings = Settings()
        dftb_settings.input = settings.specific.dftb
        result = plams.DFTBJob(name=job_name, molecule=mol,
                               settings=dftb_settings).run()

        return DFTB_Result(dftb_settings, mol, result.job.name,
                           plams_dir=result.job.path)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        msg = 'Generic keyword "' + key + '" not implemented for package DFTB.'
        warn(msg)


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results"""

    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 project_name=None):
        properties = 'data/dictionaries/propertiesDFTB.json'
        super().__init__(settings, molecule, job_name, plams_dir=plams_dir,
                         properties=properties)
        kf_filename = join(plams_dir, '{}.rkf'.format(job_name))
        self.kf = plams.kftools.KFFile(kf_filename)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        return DFTB_Result(settings, molecule, job_name,
                           archive["plams_dir"].path, project_name)

    @property
    def molecule(self, unit='bohr', internal=False, n=1):
        """WARNING: Cheap copy from PLAMS, do not keep this!!!"""
        m = self._molecule.copy()
        natoms = len(m)
        coords = self.kf.read('Molecule', 'Coords')
        coords = [coords[i: i + 3] for i in range(0, len(coords), 3)]

        if len(coords) > natoms:
            coords = coords[(n - 1) * natoms: n * natoms]
        if internal:
            mapping = self._int2inp()
            coords = [coords[mapping[i] - 1] for i in range(len(coords))]
        for at, coord in zip(m, coords):
            at.move_to(coord, 'bohr')

        return m

adf = ADF()
dftb = DFTB()
