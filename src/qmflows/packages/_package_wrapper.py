"""A module which adds the :class:`PackageWrapper` class.

The herein implemented class serves as a wrapper around the qmflows
:class:`~qmflows.packages.Package` class,
taking a single :class:`plams.Job<scm.plams.core.basejob.Job>` type as argument
and, upon calling :class:`pkg = PackageWrapper(...); pkg()<PackageWrapper.__call__>` an
appropiate instance of Package subclas instance is called.

For example, passing :class:`plams.ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>` will
automatically call :data:`~qmflows.adf`,
:class:`plams.Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>` will
call :data:`~qmflows.cp2k`, *etc*.

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

.. testsetup:: python

    >>> from scm.plams import from_smiles

    >>> mol = from_smiles('C')
    >>> settings = Settings()

.. code:: python

    >>> from scm.plams import AMSJob, Molecule, Settings

    >>> from qmflows import run, PackageWrapper
    >>> from qmflows.packages import ResultWrapper

    >>> mol = Molecule(...)  # doctest: +SKIP
    >>> settings = Settings(...)  # doctest: +SKIP

    >>> pkg = PackageWrapper(AMSJob)
    >>> job = pkg(settings, mol, name='amsjob')
    >>> result: ResultWrapper = run(job)

    >>> energy = result.get_energy()  # doctest: +SKIP
    >>> mol = result.get_molecule()  # doctest: +SKIP
    >>> freq = result.get_frequencies()  # doctest: +SKIP


Index
-----
.. currentmodule:: qmflows.packages
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
.. data:: JOB_MAP
    :annotation: : dict[type[plams.Job], Package]

    A dictionary mapping PLAMS :class:`Job<scm.plams.core.basejob.Job>` types
    to appropiate QMFlows :class:`~qmflows.packages.Package` instance.

    .. code:: python

        >>> from __future__ import annotations

        >>> from qmflows.packages import Package, cp2k, adf, dftb, orca
        >>> from scm import plams

        >>> plams.Job = plams.core.basejob.Job
        >>> plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

        >>> JOB_MAP: dict[type[plams.Job], Package] = {
        ...     plams.Cp2kJob: cp2k,
        ...     plams.ADFJob: adf,
        ...     plams.DFTBJob: dftb,
        ...     plams.ORCAJob: orca
        ... }

"""

from __future__ import annotations

import os
from os.path import join
from typing import (
    TypeVar, ClassVar, Any, Generic, TYPE_CHECKING
)
from warnings import warn

from scm import plams
from noodles import has_scheduled_methods, schedule

from ._packages import Package, Result, load_properties
from ._scm import adf, dftb
from ._orca import orca
from ._cp2k import cp2k
from .._settings import Settings
from ..type_hints import _Settings, MolType
from ..warnings_qmflows import Key_Warning

if TYPE_CHECKING:
    PT = TypeVar('PT', bound="PackageWrapper")

plams.Job = plams.core.basejob.Job
plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

__all__ = ['PackageWrapper', 'ResultWrapper', 'JOB_MAP']

JT = TypeVar("JT", bound=plams.core.basejob.Job)

JOB_MAP: dict[type[plams.Job], Package] = {
    plams.Cp2kJob: cp2k,
    plams.ADFJob: adf,
    plams.DFTBJob: dftb,
    plams.ORCAJob: orca
}


class ResultWrapper(Result):
    """The matching :class:`~qmflows.packages.Result` subclass for :class:`PackageWrapper`."""  # noqa

    prop_mapping: ClassVar[_Settings] = load_properties('PackageWrapper', prefix='properties')


@has_scheduled_methods
class PackageWrapper(Package, Generic[JT]):
    """A :class:`~qmflows.packages.Package` subclass for processing arbitrary :class:`plams.Job<scm.plams.core.basejob.Job>` types.

    Will automatically convert the passed Job type into the appropiate
    Package instance upon calling :meth:`PackageWrapper.__call__`.

    Examples
    --------
    .. code:: python

        >>> from scm.plams import ADFJob, AMSJob

        >>> from qmflows import PackageWrapper, run
        >>> from qmflows.packages import ResultWrapper, ADF_Result

        # Start of with two PackageWrapper instances
        >>> pkg_adf = PackageWrapper(ADFJob)
        >>> pkg_ams = PackageWrapper(AMSJob)

        # End up with two different Result instances
        >>> result_adf: ADF_Result = run(pkg_adf(...), ...)  # doctest: +SKIP
        >>> result_ams: ResultWrapper = run(pkg_ams(...), ...)  # doctest: +SKIP

    Attributes
    ----------
    job_type : :class:`type[plams.Job] <type>`
        The to-be executed :class:`plams.Job<scm.plams.core.basejob.Job>` type.
        Will be automatically translated, when possible, into to the appropiate
        :class:`~qmflows.packages.Package` instance upon calling :meth:`PackageWrapper.__call__`.
        If not, default to the more bare-bones implementation within this class
        and the matching :class:`ResultWrapper` instance.

    See Also
    --------
    :data:`~qmflows.packages.JOB_MAP`
        A dictionary mapping PLAMS Job types to appropiate QMFlows Package instances.

    """  # noqa

    generic_mapping: ClassVar[_Settings] = load_properties('PackageWrapper', prefix='generic2')
    result_type: ClassVar[type[ResultWrapper]] = ResultWrapper
    job_type: type[JT]

    if TYPE_CHECKING:
        def __getattr__(self, name: str) -> Any: ...

    def __init__(self, job_type: type[JT], name: None | str = None) -> None:
        """Initialize this instance.

        Parameters
        ----------
        job_type : :class:`type[plams.Job] <type>`
            The to-be executed :class:`plams.Job<scm.plams.core.basejob.Job>` type.
            Will be automatically translated, when possible, into to the appropiate
            :class:`~qmflows.packages.Package` instance upon
            calling :meth:`PackageWrapper.__call__`.
            If not, default to the more bare-bones implementation within this class
            and the matching :class:`ResultWrapper` instance.
            See also :attr:`PackageWrapper.job_type`.

        """
        if name is None:
            pkg_name = job_type.__class__.__name__.lower().rstrip('job')
        else:
            pkg_name = name
        super().__init__(pkg_name)
        self.job_type = job_type

    def __reduce__(self: PT) -> tuple[type[PT], tuple[type[JT], str]]:
        """A helper function for :mod:`pickle`."""
        return type(self), (self.job_type, self.pkg_name)

    @schedule(display="Running {self.pkg_name} {job_name}...", store=True, confirm=True)
    def __call__(self, settings: Settings,
                 mol: MolType,
                 job_name: str = '', **kwargs: Any) -> Result:
        """If possible, call :meth:`__call__()` of the Package instance appropiate to :attr:`PackageWrapper.job_type`.

        If not, default to the base :meth:`~qmflows.packages.Package.__call__` method.

        """  # noqa
        try:  # Call one of the dedicated Package subclasses
            package_obj = JOB_MAP[self.job_type]
            return package_obj(settings, mol, job_name, **kwargs)
        except KeyError:  # Continue with this instance's __call__
            return super().__call__(settings, mol, job_name, **kwargs)

    def run_job(self, settings: Settings, mol: plams.Molecule,
                job_name: str = 'job',
                work_dir: None | str | os.PathLike[str] = None,
                validate_output: bool = True,
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

        return self.result_type(
            settings, mol, job_name,
            dill_path=dill_path, plams_dir=r.job.path,
            work_dir=work_dir, status=job.status, warnings=None
        )

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Method not implemented."""
        warn('No generic keywords implemented for PackageWrapper', category=Key_Warning)
