"""A set of modules for managing various quantum-chemical packages."""

from typing import TYPE_CHECKING

from ._packages import (
    Package, Result, run, registry,
)

from ._cp2k import CP2K_Result, CP2K, cp2k
from ._cp2k_mm import CP2KMM_Result, CP2KMM, cp2k_mm
from ._scm import ADF_Result, DFTB_Result, ADF, DFTB, adf, dftb
from ._orca import ORCA_Result, ORCA, orca
from ._package_wrapper import PackageWrapper, ResultWrapper, JOB_MAP

__all__ = [
    'Package', 'Result', 'run', 'registry',
    'CP2K_Result', 'CP2K', 'cp2k',
    'CP2KMM_Result', 'CP2KMM', 'cp2k_mm',
    'ADF_Result', 'DFTB_Result', 'ADF', 'DFTB', 'adf', 'dftb',
    'ORCA_Result', 'ORCA', 'orca',
    'PackageWrapper', 'ResultWrapper', 'JOB_MAP',
]

if not TYPE_CHECKING:
    from scm import plams

    #: Placeholder docstring for sphinx.
    JOB_MAP: "dict[type[plams.Job], Package]"

del TYPE_CHECKING
