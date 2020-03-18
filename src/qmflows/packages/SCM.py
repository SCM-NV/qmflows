__all__ = ['ADF_Result', 'DFTB_Result', 'adf', 'dftb']

# =======>  Standard and third party Python Libraries <======

import os
import struct
from os.path import join
from warnings import warn
from typing import Optional, Union, Any, ClassVar

from scm import plams

from .packages import Package, package_properties, Result, get_tmpfile_name
from ..settings import Settings
from ..type_hints import WarnMap, Final
from ..warnings_qmflows import Key_Warning, QMFlows_Warning

# ========================= ADF ============================


class ADF(Package):
    """:class:`~qmflows.packages.packages.Package` subclass for ADF.

    This class takes care of calling the *ADF* quantum package.
    it uses both Plams and the Templates module to create the input
    files, while Plams takes care of running
    the :class:`plams.ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>`.
    It returns a :class:`ADF_Result` instance containing the output data.

    """

    generic_dict_file: ClassVar[str] = 'generic2ADF.yaml'

    def __init__(self) -> None:
        super(ADF, self).__init__("adf")

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_name: str = 'ADFjob', nproc: Optional[int] = None,
                **kwargs: Any) -> 'ADF_Result':
        """Execute ADF job.

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
        relative_path_t21 = join(*str(path_t21).split(os.sep)[-3:])

        # Relative job path
        relative_plams_path = join(*str(result.job.path).split(os.sep)[-2:])

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        adf_result = ADF_Result(
            adf_settings, mol, result.job.name, relative_path_t21, dill_path,
            plams_dir=relative_plams_path, status=job.status)

        return adf_result

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Handle special ADF-specific keywords.

        Some keywords provided by the user do not have a straightforward
        translation to *ADF* input and require some hooks that handles the
        special behaviour of the following keywords:

        * ``freeze``
        * ``selected_atoms``

        """
        def freeze() -> None:
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

        def selected_atoms() -> None:
            settings.specific.adf.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = f'selected_atoms {value} is not a list'
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

        def inithess() -> None:
            hess_path = get_tmpfile_name()
            hess_file = open(hess_path, "w")
            hess_file.write(" ".join(['{:.6f}'.format(v) for v in value]))
            settings.specific.adf.geometry.inithess = hess_path

        def constraint() -> None:
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
                        warn(f'Invalid constraint key: {k}', category=Key_Warning)

        # Available translations
        functions = {'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'inithess': inithess,
                     'constraint': constraint}
        if key in functions:
            functions[key]()
        else:
            warn(f'Generic keyword {key!r} not implemented for package ADF',
                 category=Key_Warning)


class ADF_Result(Result):
    """Class providing access to PLAMS ADFJob result results."""

    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 path_t21: Union[str, os.PathLike],
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 status: str = 'done',
                 warnings: Optional[WarnMap] = None) -> None:
        # Load available property parser from yaml file.
        super().__init__(settings, molecule, job_name, dill_path,
                         plams_dir=plams_dir, properties=package_properties['adf'],
                         status=status, warnings=warnings)

        # Create a KF reader instance
        self.kf = plams.KFFile(path_t21)

    def get_property_kf(self, prop: str, section: Optional[str] = None) -> Any:
        """Interface for :meth:`plams.KFFile.read()<scm.plams.tools.kftools.KFFile.read>`."""
        return self.kf.read(section, prop)

    @property
    def molecule(self, unit: str = 'bohr',
                 internal: bool = False,
                 n: int = 1) -> Optional[plams.Molecule]:
        """WARNING: Cheap copy from PLAMS, do not keep this!!!."""
        try:
            m = self._molecule.copy()
        except AttributeError:
            return None  # self._molecule can be None
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
    """:class:`~qmflows.packages.packages.Package` subclass for DFTB."""

    generic_dict_file: ClassVar[str] = 'generic2DFTB.yaml'

    def __init__(self) -> None:
        super().__init__("dftb")

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_name: str = 'DFTBjob',
                nproc: Optional[int] = None,
                **kwargs: Any) -> 'DFTB_Result':
        """Execute an DFTB job with the *ADF* quantum package.

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
        job = plams.DFTBJob(name=job_name, molecule=mol, settings=dftb_settings)

        # Check RKF status
        try:
            result = job.run()
            name = result.job.name
            path = result.job.path
        except struct.error:
            job.status = 'failed'
            name = job_name
            path = None
            warn(f"job:{job_name} has failed.\nRKF is corrupted", category=QMFlows_Warning)

        if job.status in ['failed', 'crashed']:
            plams.config.default_jobmanager.remove_job(job)

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        return DFTB_Result(dftb_settings, mol, name, dill_path,
                           plams_dir=path, status=job.status)

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Translate generic keywords to their corresponding Orca keywords."""
        def freeze() -> None:
            settings.specific.dftb.geometry.optim = "cartesian"
            if not isinstance(value, list):
                raise RuntimeError(f'freeze {value} is not a list')
            if isinstance(value[0], int):
                for a in value:
                    at = 'atom ' + str(a)
                    settings.specific.dftb.constraints[at] = ""
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol in value:
                        at = 'atom ' + str(a + 1)
                        settings.specific.dftb.constraints[at] = ""

        def selected_atoms() -> None:
            settings.specific.dftb.geometry.optim = "cartesian"
            if not isinstance(value, list):
                raise RuntimeError(f'selected_atoms {value} is not a list')
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

        def constraint() -> None:
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
                        warn(f'Invalid constraint key: {k}', category=Key_Warning)

        # Available translations
        functions = {'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'constraint': constraint}
        if key in functions:
            functions[key]()
        else:
            warn(f'Generic keyword {key!r} not implemented for package DFTB',
                 category=Key_Warning)


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results."""

    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 status: str = 'done',
                 warnings: Optional[WarnMap] = None) -> None:
        # Read available propiety parsers from a yaml file
        super().__init__(settings, molecule, job_name, dill_path,
                         plams_dir=plams_dir, properties=package_properties['dftb'],
                         status=status, warnings=warnings)

        if plams_dir is not None:
            kf_filename = join(plams_dir, f'{job_name}.rkf')
            # create a kf reader instance
            self.kf = plams.KFFile(kf_filename)
        else:
            self.kf = None

    @property
    def molecule(self, unit: str = 'bohr',
                 internal: bool = False,
                 n: int = 1) -> plams.Molecule:
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


#: An instance :class:`ADF`.
#: Only one instance of this class should exist at any given momemt;
#: *i.e.* this value is a singleton.
adf: Final[ADF] = ADF()

#: An instance :class:`DFTB`.
#: Only one instance of this class should exist at any given momemt;
#: *i.e.* this value is a singleton.
dftb: Final[DFTB] = DFTB()
