__all__ = ['adf', 'dftb']

# =======>  Standard and third party Python Libraries <======

from os.path import join
from warnings import warn
from qmflows.settings import Settings
from qmflows.packages.packages import (Package, package_properties, Result, get_tmpfile_name)
from scm import plams

import builtins
import struct

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
    def run_job(settings, mol, job_name='ADFjob', nproc=None):
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
        :returns: :class:`~qmflows.packages.SCM.ADF_Result`
        """
        adf_settings = Settings()
        if nproc:
            adf_settings.runscript.nproc = nproc
        adf_settings.input = settings.specific.adf
        job = plams.ADFJob(name=job_name, molecule=mol,
                                               settings=adf_settings)
        result = job.run()
        # Path to the tape 21 file
        path_t21 = result._kf.path

        # Relative path to the CWD
        relative_path_t21 = '/'.join(path_t21.split('/')[-3:])

        # Relative job path
        relative_plams_path = '/'.join(result.job.path.split('/')[-2:])

        adf_result = ADF_Result(
            adf_settings, mol, result.job.name, relative_path_t21,
            plams_dir=relative_plams_path, status=job.status)

        return adf_result

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
        def freeze():
            settings.specific.adf.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'freeze ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in value:
                    at = 'atom ' + str(a)
                    settings.specific.adf.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[at] = ""

        def selected_atoms():
            settings.specific.adf.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in range(len(mol)):
                    if a + 1 not in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol not in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.adf.constraints[at] = ""

        def inithess():
            hess_path = get_tmpfile_name()
            hess_file = open(hess_path, "w")
            hess_file.write(" ".join(['{:.6f}'.format(v) for v in value]))
            settings.specific.adf.geometry.inithess = hess_path

        def constraint():
            if isinstance(value, Settings):
                for k, v in value.items():
                    ks = k.split()
                    # print('--->', ks, type(ks[2]), type(value), v)
                    if ks[0] == 'dist' and len(ks) == 3:
                        name = 'dist {:d} {:d}'.format(int(ks[1]),
                                                       int(ks[2]))
                        settings.specific.adf.constraints[name] = v
                    elif ks[0] == 'angle' and len(ks) == 4:
                        name = 'angle {:d} {:d} {:d}'.format(int(ks[1]),
                                                             int(ks[2]),
                                                             int(ks[3]))
                        settings.specific.adf.constraints[name] = v
                    elif ks[0] == 'dihed' and len(ks) == 5:
                        name = 'dihed {:d} {:d} {:d} {:d}'.\
                            format(int(ks[1]), int(ks[2]),
                                   int(ks[3]), int(ks[4]))
                        settings.specific.adf.constraints[name] = v
                    else:
                        warn('Invalid constraint key: ' + k)

        # Available translations
        functions = {'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'inithess': inithess,
                     'constraint': constraint}
        if key in functions:
            functions[key]()
        else:
            msg = 'Generic keyword "' + key + '" not implemented for package ADF.'
            warn(msg)


class ADF_Result(Result):
    """Class providing access to PLAMS ADFJob result results"""

    def __init__(self, settings, molecule, job_name, path_t21, plams_dir=None,
                 status='done', warnings=None):
        # Load available property parser from Json file.
        properties = package_properties['adf']
        super().__init__(settings, molecule, job_name, plams_dir=plams_dir,
                         properties=properties, status=status, warnings=warnings)
        # Create a KF reader instance
        self.kf = plams.KFFile(path_t21)

    def __deepcopy__(self, memo):
        return ADF_Result(self.settings,
                          self._molecule,
                          self.job_name,
                          self.kf.path,
                          plams_dir = self.archive['plams_dir'].path,
                          status = self.status,
                          warnings = self.warnings
                          )

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status, warnings):
        """
        Methods to deserialize an `ADF_Result` object.
        """
        plams_dir = archive["plams_dir"].path
        path_t21 = join(plams_dir, '{}.t21'.format(job_name))
        return ADF_Result(settings, molecule, job_name, path_t21, plams_dir,
                          status, warnings)

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
    def run_job(settings, mol, job_name='DFTBjob', nproc=None):
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
        :returns: :class:`~qmflows.packages.SCM.DFTB_Result`
        """
        dftb_settings = Settings()
        if nproc:
            dftb_settings.runscript.nproc = nproc
        dftb_settings.input = settings.specific.dftb
        job = plams.DFTBJob(name=job_name, molecule=mol,
                                                settings=dftb_settings)

        # Check RKF status
        try:
            result = job.run()
            name = result.job.name
            path = result.job.path
        except struct.error:
            job.status = 'failed'
            name = job_name
            path = None
            msg = "job:{} has failed.\nRKF is corrupted"
            print(msg.format(job_name))

        if job.status in ['failed', 'crashed']:
            builtins.config.jm.remove_job(job)

        return DFTB_Result(dftb_settings, mol, name,
                           plams_dir=path, status=job.status)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Translate generic keywords to their corresponding Orca keywords.
        """
        def freeze():
            settings.specific.dftb.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'freeze ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in value:
                    at = 'atom ' + str(a)
                    settings.specific.dftb.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.dftb.constraints[at] = ""

        def selected_atoms():
            settings.specific.dftb.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                for a in range(len(mol)):
                    if a + 1 not in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.dftb.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol not in value:
                        name = 'atom ' + str(a + 1)
                        settings.specific.dftb.constraints[name] = ""

        def constraint():
            if isinstance(value, Settings):
                for k, v in value.items():
                    ks = k.split()
                    if ks[0] == 'dist' and len(ks) == 3:
                        name = 'dist {:d} {:d}'.format(int(ks[1]),
                                                       int(ks[2]))
                        settings.specific.dftb.constraints[name] = v
                    elif ks[0] == 'angle' and len(ks) == 4:
                        name = 'angle {:d} {:d} {:d}'.format(
                            int(ks[1]), int(ks[2]), int(ks[2]))
                        settings.specific.dftb.constraints[name] = v
                    elif ks[0] == 'dihed' and len(ks) == 5:
                        name = 'dihed {:d} {:d} {:d} {:d}'.format(
                            int(ks[1]), int(ks[2]),
                            int(ks[3]), int(ks[4]))
                        settings.specific.dftb.constraints[name] = v
                    else:
                        warn('Invalid constraint key: ' + k)

        # Available translations
        functions = {'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'constraint': constraint}
        if key in functions:
            functions[key]()
        else:
            msg = 'Generic keyword "' + key + '" not implemented for package DFTB.'
            warn(msg)


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results"""

    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 status='done', warnings=None):
        # Read available propiety parsers from a JSON file
        properties = package_properties['dftb']
        super().__init__(settings, molecule, job_name, plams_dir=plams_dir,
                         properties=properties, status=status, warnings=warnings)
        if plams_dir is not None:
            kf_filename = join(plams_dir, '{}.rkf'.format(job_name))
            # create a kf reader instance
            self.kf = plams.KFFile(kf_filename)
        else:
            self.kf = None

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status, warnings):
        return DFTB_Result(settings, molecule, job_name,
                           archive["plams_dir"].path, status, warnings)

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
