"""Utility functions for the QMFlows tests.

Index
-----
.. currentmodule:: qmflows.test_utils
.. autosummary::
    fill_cp2k_defaults
    get_mm_settings
    PATH
    PATH_MOLECULES
    HAS_RDKIT
    requires_cp2k
    requires_orca
    requires_ams
    requires_adf

API
---
.. autofunction:: fill_cp2k_defaults
.. autofunction:: get_mm_settings
.. autodata:: PATH
    :annotation: : pathlib.Path
.. autodata:: PATH_MOLECULES
    :annotation: : pathlib.Path
.. autodata:: HAS_RDKIT
    :annotation: : bool
.. autodata:: requires_cp2k
.. autodata:: requires_orca
.. autodata:: requires_ams
.. autodata:: requires_adf

"""

import os
from pathlib import Path

import pytest
from distutils.spawn import find_executable

from .fileFunctions import yaml2Settings
from .settings import Settings
from .warnings_qmflows import Assertion_Warning

__all__ = [
    'get_mm_settings',
    'PATH',
    'PATH_MOLECULES',
    'Assertion_Warning',
    'HAS_RDKIT',
    'requires_cp2k',
    'requires_orca'
    'requires_ams',
    'requires_adf',
]

try:
    import rdkit
except ImportError:
    #: Whether RDKit has been installed.
    HAS_RDKIT = False
else:
    del rdkit
    HAS_RDKIT = True

#: The path to the ``tests/test_files`` directory.
PATH = Path('test') / 'test_files'

#: The path to the ``tests/test_files/molecules`` directory.
PATH_MOLECULES = PATH / "molecules"

#: Basisset specification for CP2K calculations
kinds_template = yaml2Settings("""
specific:
  cp2k:
     force_eval:
       subsys:
         kind:
           C:
             basis_set: DZVP-MOLOPT-SR-GTH-q4
             potential: GTH-PBE-q4
           H:
             basis_set: DZVP-MOLOPT-SR-GTH-q1
             potential: GTH-PBE-q1
           O:
             basis_set: DZVP-MOLOPT-SR-GTH-q6
             potential: GTH-PBE-q6
""")


def fill_cp2k_defaults(s: Settings) -> Settings:
    """Fill missing values from a job template."""
    s.periodic = "None"
    s.cell_parameters = 10
    s = s.overlay(kinds_template)

    # functional
    s.specific.cp2k.force_eval.dft.xc.xc_functional.pbe = {}

    # basis and potential
    add_basis_potential(s)

    return s


def add_basis_potential(s: Settings) -> None:
    """Add basis and potential path to the settings."""
    s.specific.cp2k.force_eval.dft.potential_file_name = (
        PATH / "GTH_POTENTIALS").absolute().as_posix()
    s.specific.cp2k.force_eval.dft.basis_set_file_name = (
        PATH / "BASIS_MOLOPT").absolute().as_posix()


def get_mm_settings() -> Settings:
    """Construct and return CP2K settings for classical forcefield calculations.

    Note that these settings will still have to be updated with a job-specific template.

    """
    charge = Settings()
    charge.param = 'charge'
    charge.Cd = 0.9768
    charge.Se = -0.9768
    charge.O2D2 = -0.4704
    charge.C2O3 = 0.4524

    lj = Settings()
    lj.param = 'epsilon', 'sigma'
    lj.unit = 'kcalmol', 'angstrom'

    lj['Cd Cd  '] = 0.0741, 1.2340
    lj['Cd O2D2'] = 0.4383, 2.4710
    lj['Cd Se  '] = 0.3639, 2.9400
    lj['Se O2D2'] = 0.3856, 3.5260
    lj['Se Se  '] = 0.1020, 4.8520

    lj['Cd C331'] = lj['Cd C2O3'] = 0.1547, 2.9841
    lj['Se C331'] = lj['Se C2O3'] = 0.1748, 3.5885
    lj['Cd HGA3'] = 0.1002, 2.5542
    lj['Se HGA3'] = 0.1132, 3.1587

    s = Settings()
    s.psf = PATH / 'Cd68Cl26Se55__26_acetate.psf'
    s.prm = PATH / 'Cd68Cl26Se55__26_acetate.prm'
    s.charge = charge
    s.lennard_jones = lj
    s.periodic = 'none'
    s.cell_parameters = [50, 50, 50]
    return s


def _has_exec(executable: str) -> bool:
    """Check if the passed executable is installed."""
    path = find_executable(executable)
    return path is not None


def _has_env_vars(*env_vars: str) -> bool:
    """Check if the passed environment variables are available."""
    return set(env_vars).issubset(os.environ)


#: A mark for skipping tests if CP2K is not installed
requires_cp2k = pytest.mark.skipif(not _has_exec("cp2k.popt"), reason="Requires CP2K")

#: A mark for skipping tests if Orca is not installed
requires_orca = pytest.mark.skipif(not _has_exec("orca"), reason="Requires Orca")

#: A mark for skipping tests if AMS >=2020 is not installed
requires_ams = pytest.mark.skipif(
    not _has_env_vars("AMSBIN", "AMSHOME", "AMSRESOURCES"),
    reason="Requires AMS >=2020",
)

#: A mark for skipping tests if ADF <=2019 is not installed
requires_adf = pytest.mark.skipif(
    not _has_env_vars("ADFBIN", "ADFHOME", "ADFRESOURCES"),
    reason="Requires ADF <=2019",
)
