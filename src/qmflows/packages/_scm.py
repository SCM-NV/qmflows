"""Interface to the SCM funcionality."""

from __future__ import annotations

import os
import struct
from os.path import join, basename, normpath
from typing import Any, ClassVar, Final
from warnings import warn

import numpy as np
from scm import plams

from ._packages import Package, Result, load_properties
from .._settings import Settings
from ..type_hints import WarnMap, _Settings
from ..warnings_qmflows import Key_Warning, QMFlows_Warning
from ..utils import get_tmpfile_name

__all__ = ['ADF_Result', 'DFTB_Result', 'ADF', 'DFTB', 'adf', 'dftb']

# ========================= ADF ============================


def handle_SCM_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule, package: str) -> None:
    """Handle special ADF/DFTB specific keywords."""
    def freeze() -> None:
        settings.specific[package].geometry.optim = "cartesian"
        if not isinstance(value, list):
            msg = 'freeze ' + str(value) + ' is not a list'
            raise RuntimeError(msg)
        if isinstance(value[0], int):
            for a in value:
                at = 'atom ' + str(a)
                settings.specific[package].constraints[at] = ""
        else:
            for a in range(len(mol)):
                if mol[a + 1].symbol in value:
                    at = 'atom ' + str(a + 1)
                    settings.specific[package].constraints[at] = ""

    def selected_atoms() -> None:
        settings.specific[package].geometry.optim = "cartesian"
        if not isinstance(value, list):
            msg = f'selected_atoms {value} is not a list'
            raise RuntimeError(msg)
        if isinstance(value[0], int):
            for a in range(len(mol)):
                if a + 1 not in value:
                    at = 'atom ' + str(a + 1)
                    settings.specific[package].constraints[at] = ""
        else:
            for a in range(len(mol)):
                if mol[a + 1].symbol not in value:
                    at = 'atom ' + str(a + 1)
                    settings.specific[package].constraints[at] = ""

    def inithess() -> None:
        hess_path = get_tmpfile_name("ADF_hessian_")
        hess_file = open(hess_path, "w")
        hess_file.write(" ".join(['{:.6f}'.format(v) for v in value]))
        settings.specific[package].geometry.inithess = hess_path.as_posix()

    def constraint() -> None:
        if isinstance(value, Settings):
            for k, v in value.items():
                ks = k.split()
                # print('--->', ks, type(ks[2]), type(value), v)
                if ks[0] == 'dist' and len(ks) == 3:
                    name = 'dist {:d} {:d}'.format(int(ks[1]),
                                                   int(ks[2]))
                    settings.specific[package].constraints[name] = v
                elif ks[0] == 'angle' and len(ks) == 4:
                    name = 'angle {:d} {:d} {:d}'.format(int(ks[1]),
                                                         int(ks[2]),
                                                         int(ks[3]))
                    settings.specific[package].constraints[name] = v
                elif ks[0] == 'dihed' and len(ks) == 5:
                    name = 'dihed {:d} {:d} {:d} {:d}'.\
                        format(int(ks[1]), int(ks[2]),
                               int(ks[3]), int(ks[4]))
                    settings.specific[package].constraints[name] = v
                else:
                    warn(
                        f'Invalid constraint key: {k}', category=Key_Warning)

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

    prop_mapping: ClassVar[_Settings] = load_properties('ADF', prefix='properties')

    # Attributes accessed via `__getattr__`
    charges: None | Any
    dipole: None | Any
    energy: None | Any
    enthalpy: None | Any
    free_energy: None | Any
    frequencies: None | Any
    hessian: None | Any
    homo: None | Any
    lumo: None | Any
    optcycles: None | Any
    runtime: None | Any

    def __init__(self, settings: None | Settings,
                 molecule: None | plams.Molecule,
                 job_name: str,
                 dill_path: None | str | os.PathLike[str] = None,
                 plams_dir: None | str | os.PathLike[str] = None,
                 work_dir: None | str | os.PathLike[str] = None,
                 status: str = 'successful',
                 warnings: None | WarnMap = None) -> None:
        # Load available property parser from yaml file.
        super().__init__(settings, molecule, job_name, dill_path,
                         plams_dir=plams_dir, status=status, warnings=warnings)

        # Create a KF reader instance
        if work_dir is not None and plams_dir is not None:
            # The t21 path has to be absolute: use workdir instead of plams_dir
            name_t21 = basename(normpath(work_dir))
            path_t21 = join(plams_dir, f'{name_t21}.t21')
            self.kf: "None | plams.KFFile" = plams.KFFile(path_t21)
        else:
            self.kf = None

    def get_property_kf(self, prop: str, section: None | str = None) -> None | Any:
        """Interface for :meth:`plams.KFFile.read()<scm.plams.tools.kftools.KFFile.read>`."""
        if self.kf is None:
            return None
        else:
            return self.kf.read(section, prop)

    @property
    def molecule(self) -> None | plams.Molecule:
        """WARNING: Cheap copy from PLAMS, do not keep this!!!."""
        if self._molecule is None or self.kf is None:
            return None
        else:
            m = self._molecule.copy()
        # Find out correct location
        coords = self.kf.read('Geometry', 'xyz InputOrder')
        coords = np.array([coords[i: i + 3] for i in range(0, len(coords), 3)])
        m.from_array(coords)
        return m

    @property
    def geometry(self) -> None | plams.Molecule:
        """An alias for :attr:`ADF_Result.molecule`."""
        return self.molecule


class DFTB_Result(Result):
    """Class providing access to PLAMS DFTBJob result results."""

    prop_mapping: ClassVar[_Settings] = load_properties('DFTB', prefix='properties')

    # Attributes accessed via `__getattr__`
    charges: "None | Any"
    dipole: "None | Any"
    energy: "None | Any"
    enthalpy: "None | Any"
    free_energy: "None | Any"
    frequencies: "None | Any"
    hessian: "None | Any"

    def __init__(self, settings: None | Settings,
                 molecule: None | plams.Molecule,
                 job_name: str,
                 dill_path: None | str | os.PathLike[str] = None,
                 plams_dir: None | str | os.PathLike[str] = None,
                 work_dir: None | str | os.PathLike[str] = None,
                 status: str = 'successful',
                 warnings: None | WarnMap = None) -> None:
        # Read available propiety parsers from a yaml file
        super().__init__(settings, molecule, job_name, dill_path,
                         plams_dir=plams_dir, status=status, warnings=warnings)

        if plams_dir is not None:
            kf_filename = join(plams_dir, 'dftb.rkf')
            # create a kf reader instance
            self.kf: "None | plams.KFFile" = plams.KFFile(kf_filename)
        else:
            self.kf = None

    @property
    def molecule(self) -> None | plams.Molecule:
        """Read molecule from output."""
        if self._molecule is None or self.kf is None:
            return None
        else:
            m = self._molecule.copy()
        coords = self.kf.read('Molecule', 'Coords')
        coords = np.array([coords[i: i + 3] for i in range(0, len(coords), 3)])
        m.from_array(coords)
        return m

    @property
    def geometry(self) -> None | plams.Molecule:
        """An alias for :attr:`DFTB_Result.molecule`."""
        return self.molecule


class ADF(Package):
    """:class:`~qmflows.packages.Package` subclass for ADF.

    This class takes care of calling the *ADF* quantum package.
    it uses both Plams and the Templates module to create the input
    files, while Plams takes care of running
    the :class:`plams.ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>`.
    It returns a :class:`ADF_Result` instance containing the output data.

    """

    generic_mapping: ClassVar[_Settings] = load_properties('ADF', prefix='generic2')
    result_type: ClassVar[type[ADF_Result]] = ADF_Result

    def __init__(self, pkg_name: str = "adf") -> None:
        super().__init__(pkg_name)

    @classmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule,
                job_name: str = 'ADFjob', nproc: None | int = None,
                validate_output: bool = True,
                **kwargs: Any) -> ADF_Result:
        """Execute ADF job.

        :param settings: user input settings.
        :type settings: :class:`~qmflows.Settings`
        :param mol: Molecule to run the simulation
        :type mol: Plams Molecule
        :parameter input_file_name: The user can provide a name for the
                                    job input.
        :type input_file_name: String
        :parameter out_file_name: The user can provide a name for the
                                  job output.
        :type out_file_name: String
        :returns: :class:`~qmflows.packages.ADF_Result`

        """
        adf_settings = Settings()
        if nproc:
            adf_settings.runscript.nproc = nproc
        adf_settings.input = settings.specific.adf
        job = plams.ADFJob(name=job_name, molecule=mol,
                           settings=adf_settings)
        result = job.run()

        # Relative job path
        relative_plams_path = join(*str(result.job.path).split(os.sep)[-2:])

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        adf_result = cls.result_type(
            adf_settings, mol, result.job.name, dill_path,
            plams_dir=relative_plams_path,
            work_dir=result.job.path, status=job.status)

        return adf_result

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Handle special ADF/DFTB specific keywords.

        Some keywords provided by the user do not have a straightforward
        translation to *ADF* input and require some hooks that handles the
        special behaviour of the following keywords:

        * ``freeze``
        * ``selected_atoms``
        * ``initHess``
        * ``Constraint``
        """
        handle_SCM_special_keywords(settings, key, value, mol, "adf")


class DFTB(Package):
    """:class:`~qmflows.packages.Package` subclass for DFTB."""

    generic_mapping: ClassVar[_Settings] = load_properties('DFTB', prefix='generic2')
    result_type: ClassVar[type[DFTB_Result]] = DFTB_Result

    def __init__(self, pkg_name: str = "dftb") -> None:
        super().__init__(pkg_name)

    @classmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule, job_name: str = "DFTBJob",
                work_dir: None | str | os.PathLike[str] = None,
                validate_output: bool = True,
                **kwargs: Any) -> DFTB_Result:
        """Execute an DFTB job with the AMS driver.

        In order to run a DFTB calculation we need both an AMS and DFTB sections.
        The AMS sections, specifies which tasks to run: single point, optimization, etc.
        While the DFTB section only set the input for the method.
        For more information, see: `AMS <https://www.scm.com/doc/plams/interfaces/ams.html>`
        """
        dftb_settings = Settings()
        dftb_settings.input = settings.specific.dftb
        dftb_settings.input += settings.specific.ams

        job = plams.AMSJob(name=job_name, molecule=mol,
                           settings=dftb_settings)

        # Check RKF status
        try:
            result = job.run()
            name = result.job.name
            path: None | str = result.job.path
        except struct.error:
            job.status = 'failed'
            name = job_name
            path = None
            warn(f"job:{job_name} has failed.\nRKF is corrupted",
                 category=QMFlows_Warning)

        if job.status in ['failed', 'crashed']:
            plams.config.default_jobmanager.remove_job(job)

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        return cls.result_type(dftb_settings, mol, name, dill_path,
                               plams_dir=path, status=job.status)

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Translate generic keywords to their corresponding Orca keywords."""
        handle_SCM_special_keywords(settings, key, value, mol, "dftb")


adf: Final[ADF] = ADF()
dftb: Final[DFTB] = DFTB()
