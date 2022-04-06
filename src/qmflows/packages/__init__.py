"""A set of modules for managing various quantum-chemical packages."""

from ._packages import (
    Package, Result, run,
    _load_properties as load_properties,
    _registry as registry,
)
from ._serializer import (
    _SerMolecule as SerMolecule,
    _SerSettings as SerSettings,
)

from ._cp2k import CP2K_Result, CP2K, cp2k
from ._cp2k_mm import CP2KMM_Result, CP2KMM, cp2k_mm
from ._scm import ADF_Result, DFTB_Result, ADF, DFTB, adf, dftb
from ._orca import ORCA_Result, ORCA, orca
from ._package_wrapper import PackageWrapper, ResultWrapper, JOB_MAP

__all__ = [
    'Package', 'Result', 'run',
    'CP2K_Result', 'CP2K', 'cp2k',
    'CP2KMM_Result', 'CP2KMM', 'cp2k_mm',
    'ADF_Result', 'DFTB_Result', 'ADF', 'DFTB', 'adf', 'dftb',
    'ORCA_Result', 'ORCA', 'orca',
    'PackageWrapper', 'ResultWrapper', 'JOB_MAP',
]
