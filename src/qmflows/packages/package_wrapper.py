"""A module which adds the :class:`PackageWrapper` class.

The herein implemented class serves as a wrapper around the qmflows
:class:`Result<qmflows.packages.packages.Package>` class,
taking a single :class:`plams.Job<scm.plams.core.basejob.Job>` type as argument
and, upon calling :class:`pkg = PackageWrapper(...); pkg()<PackageWrapper.__call__>` an
appropiate instance of Package subclas instance is called.

For example, passing :class:`plams.ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>` will
automatically call :data:`adf<qmflows.packages.SCM.adf>`,
:class:`plams.Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>` will
call :data:`cp2k<qmflows.packages.cp2k_package.cp2k>`, *etc*.

When no appropiate Package is found, let's say after passing the :class:`MyFancyJob` type,
the PackageWrapper class will still run the job as usual and return the matching
:class:`ResultWrapper` object.

There are however three caveats:

1. No generic keywords are implemented for such jobs.
2. Specialized warning message will not be available.
3. The availability of property-extraction methods is limited.

The therein embedded results can be still extracted by calling the
:class:`plams.Results<scm.plams.core.results.Results>` methods appropiate to the passed job type,
*e.g.* :class:`plams.AMSResults<scm.plams.interfaces.adfsuite.ams.AMSResults>`.

For example:

.. code:: python

    >>> from scm.plams import AMSJob, Molecule, Settings

    >>> from qmflows import run, PackageWrapper
    >>> from qmflows.packages.package_wrapper import ResultWrapper

    >>> mol = Molecule(...)
    >>> settings = Settings(...)

    >>> pkg = PackageWrapper(AMSJob)
    >>> job = pkg(settings, mol, name='amsjob')
    >>> result: ResultWrapper = run(job, ...)

    >>> energy = result.get_energy()  # Alias for AMSResults.get_energy()
    >>> mol = result.get_molecule()  # Alias for AMSResults.get_molecule()
    >>> freq = result.get_frequencies()  # Alias for AMSResults.get_frequencies()


Index
-----
.. currentmodule:: qmflows.packages.package_wrapper
.. autosummary::
    PackageWrapper
    PackageWrapper.__init__
    PackageWrapper.__call__
    PackageWrapper.run_job
    PackageWrapper.handle_special_keywords
    ResultWrapper
    JOB_MAP

API
---
.. autoclass:: PackageWrapper
.. automethod:: PackageWrapper.__init__
.. automethod:: PackageWrapper.__call__
.. automethod:: PackageWrapper.run_job
.. automethod:: PackageWrapper.handle_special_keywords
.. autoclass:: ResultWrapper
.. autodata:: JOB_MAP
    :annotation: : dict [type [plams.Job], Package]

    .

    .. code:: python

        >>> from typing import Dict, Type

        >>> from qmflows.packages import Package, cp2k, adf, dftb, gamess, orca
        >>> from scm import plams

        >>> plams.Job = plams.core.basejob.Job
        >>> plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

        >>> JOB_MAP: Dict[Type[plams.Job], Package] = {
        ...     plams.Cp2kJob: cp2k,
        ...     plams.ADFJob: adf,
        ...     plams.DFTBJob: dftb,
        ...     plams.GamessJob: gamess,
        ...     plams.ORCAJob: orca
        ... }

"""

import os
from os.path import join
from typing import Type, Optional, TypeVar, Union, ClassVar, Any, Dict
from warnings import warn

from scm import plams
from rdkit import Chem
from noodles import has_scheduled_methods, schedule

from .packages import Package, Result, package_properties
from .SCM import adf, dftb
from .orca import orca
from .gamess import gamess
from .cp2k_package import cp2k
from ..settings import Settings
from ..type_hints import WarnMap
from ..warnings_qmflows import Key_Warning

plams.Job = plams.core.basejob.Job
plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

__all__ = ['PackageWrapper']

#: A :class:`dict` mapping PLAMS :class:`Job<scm.plams.core.basejob.Job>` types
#: to appropiate QMFlows :class:`Package<qmflows.packages.packages.Package>` instance
JOB_MAP: Dict[Type[plams.Job], Package] = {
    plams.Cp2kJob: cp2k,
    plams.ADFJob: adf,
    plams.DFTBJob: dftb,
    plams.GamessJob: gamess,
    plams.ORCAJob: orca
}

#: TypeVar for Result objects and its subclasses.
RT = TypeVar('RT', bound=Result)


class ResultWrapper(Result):
    """The matching :class:`Result<qmflows.packages.packages.Result>` subclass for :class:`PackageWrapper`."""  # noqa

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
    """A :class:`Package<qmflows.packages.packages.Package>` subclass for processing arbitrary :class:`plams.Job<scm.plams.core.basejob.Job>` types.

    Will automatically convert the passed Job type into the appropiate
    Package instance upon calling :meth:`PackageWrapper.__call__`.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import ADFJob, AMSJob

        >>> from qmflows PackageWrapper, run
        >>> from qmflows.packages.package_wrapper import ResultWrapper
        >>> from qmflows.packages.SCM import ADF_Result

        # Start of with two PackageWrapper instances
        >>> pkg_adf = PackageWrapper(ADFJob)
        >>> pkg_ams = PackageWrapper(AMSJob)

        # End up with two different Result instances
        >>> result_adf: ADF_Result = run(pkg_adf(...), ...)
        >>> result_ams: ResultWrapper = run(pkg_ams(...), ...)

    Attributes
    ----------
    job_type : :class:`type` [:class:`plams.Job<scm.plams.core.basejob.Job>`]
        The to-be executed :class:`plams.Job<scm.plams.core.basejob.Job>` type.
        Will be automatically translated, when possible, into to the appropiate
        :class:`Package<qmflows.packages.packages.Package>` instance upon calling
        :meth:`PackageWrapper.__call__`.
        If not, default to the more bare-bones implementation within this class
        and the matching :class:`ResultWrapper` instance.

    See Also
    --------
    :data:`JOB_MAP<qmflows.packages.package_wrapper.JOB_MAP>` : :class:`dict` [:class:`type` [:class:`plams.Job<scm.plams.core.basejob.Job>`], :class:`Package<qmflows.packages.packages.Package>`]
        A :class:`dict` mapping PLAMS Job types to appropiate QMFlows Package instances.

    """  # noqa

    generic_dict_file: ClassVar[str] = 'generic2None.yaml'
    generic_package: ClassVar[bool] = True

    def __init__(self, job_type: Type[plams.Job]) -> None:
        """Initialize this instance.

        Parameters
        ----------
        job_type : :class:`type` [:class:`plams.Job<scm.plams.core.basejob.Job>`]
            The to-be executed :class:`plams.Job<scm.plams.core.basejob.Job>` type.
            Will be automatically translated, when possible, into to the appropiate
            :class:`Package<qmflows.packages.packages.Package>` instance upon calling
            :meth:`PackageWrapper.__call__`.
            If not, default to the more bare-bones implementation within this class
            and the matching :class:`ResultWrapper` instance.
            See also :attr:`PackageWrapper.job_type`.

        """
        pkg_name = job_type.__class__.__name__.lower().rstrip('job')
        super().__init__(pkg_name)
        self.job_type = job_type

    @schedule(display="Running {self.pkg_name} {job_name}...", store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: Union[plams.Molecule, Chem.Mol],
                 job_name: str = '', **kwargs: Any) -> RT:
        """If possible, call :meth:`__call__()` of the Package instance appropiate to :attr:`PackageWrapper.job_type`.

        If not, default to the base :meth:`Package.__call__()<qmflows.packages.packages.Package.__call__>` method.

        """  # noqa
        try:  # Call one of the dedicated Package subclasses
            package_obj = JOB_MAP[self.job_type]
            return package_obj(settings, mol, job_name, **kwargs)
        except KeyError:  # Continue with this instance's __call__
            return super().__call__(settings, mol, job_name, **kwargs)

    def run_job(self, settings: Settings, mol: plams.Molecule,
                job_name: str = 'job',
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> ResultWrapper:
        """Run the job and pass the resulting :class:`plams.Results<scm.plams.core.results.Results>` object to :class:`ResultWrapper`."""  # noqa
        # Input modifications
        job_settings = Settings()
        for s in settings.specific.values():
            job_settings.input.update(s)

        # Create a Plams job
        job = self.job_type(name=job_name, settings=job_settings, molecule=mol)
        r = job.run()

        # Extract the working directory
        work_dir = work_dir if work_dir is not None else job.path

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        return ResultWrapper(settings, mol, job_name, r.job.path, dill_path,
                             work_dir=work_dir, status=job.status, warnings=None)

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Method not implemented."""
        warn(f'No generic keywords implemented for PackageWrapper',
             category=Key_Warning)
