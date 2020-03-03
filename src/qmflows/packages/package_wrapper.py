import os
from os.path import join
from typing import Type, Mapping, Optional, TypeVar, Union, ClassVar, Any

from scm import plams
from rdkit import Chem
from noodles import has_scheduled_methods, schedule
from qmflows.packages.packages import Package, Result, WarnMap, package_properties
from qmflows import Settings, orca, adf, dftb, cp2k, gamess

plams.Job = plams.core.basejob.Job
plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

__all__ = []

#: Map a PLAMS Job type to an appropiate packages instance.
JOB_MAP: Mapping[Type[plams.Job], Package] = {
    plams.Cp2kJob: cp2k,
    plams.ADFJob: adf,
    plams.DFTBJob: dftb,
    plams.GamessJob: gamess,
    plams.ORCAJob: orca
}

#: TypeVar for Result objects and its subclasses.
RT = TypeVar('RT', bound=Result)


class ResultWrapper(Result):
    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 status: str = 'successful',
                 warnings: Optional[WarnMap] = None) -> None:
        """Initialize this instance."""
        super().__init__(settings, molecule, job_name, plams_dir, dill_path,
                         work_dir=work_dir, properties=package_properties[None],
                         status=status, warnings=warnings)


@has_scheduled_methods
class PackageWrapper(Package):
    generic_dict_file: ClassVar[str] = 'generic2None.json'
    generic_package: ClassVar[bool] = True

    def __init__(self, job_type: Type[plams.Job]) -> None:
        """Initialize this instance."""
        pkg_name = job_type.__class__.__name__.lower().rstrip('job')
        super().__init__(pkg_name)
        self.job_type = job_type

    @schedule(display="Running {self.pkg_name} {job_name}...",
              store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: Union[plams.Molecule, Chem.Mol],
                 job_name: str = '', **kwargs: Any) -> RT:
        """Call one of the appropiate Package.__call__ methods if possible."""
        try:
            package_obj = JOB_MAP[self.job_type]
            return package_obj(settings, mol, job_name, **kwargs)
        except KeyError:
            return super().__call__(settings, mol, job_name, **kwargs)

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_type: Type[plams.Job],
                job_name: str = 'job',
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> ResultWrapper:
        """Call the Cp2K binary using plams interface.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the numerical results.
        :type hdf5_file: String
        :param input_file_name: Optional name for the input.
        :type input_file_name: String
        :param out_file_name: Optional name for the output.
        :type out_file_name: String
        :param store_in_hdf5: wether to store the output arrays in HDF5 format.
        :type store_in_hdf5: Bool

        """
        # Input modifications
        settings = Settings()
        settings.input = settings

        # Create a Plams job
        job = job_type(name=job_name, settings=settings, molecule=mol)
        r = job.run()

        # Extract the working directory
        work_dir = work_dir if work_dir is not None else job.path

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        result = ResultWrapper(settings, mol, job_name, r.job.path, dill_path,
                               work_dir=work_dir, status=job.status, warnings=None)
        return result

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        return None
