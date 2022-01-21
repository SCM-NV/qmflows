"""Utility functions for the QMFlows tests.

Index
-----
.. currentmodule:: qmflows.test_utils
.. autosummary::
    fill_cp2k_defaults
    get_mm_settings
    validate_status
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
.. autofunction:: validate_status
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
import textwrap
from pathlib import Path

import pytest
from distutils.spawn import find_executable

from .settings import Settings
from .warnings_qmflows import Assertion_Warning
from .packages import Result

__all__ = [
    'get_mm_settings',
    'PATH',
    'PATH_MOLECULES',
    'Assertion_Warning',
    'validate_status',
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


def fill_cp2k_defaults(s: Settings) -> Settings:
    """Fill missing values from a job template.

    Returns a copy of passed template.
    """
    s = s.copy()
    s.periodic = "None"
    s.cell_parameters = 10

    # functional
    s.specific.cp2k.force_eval.dft.xc.xc_functional.pbe = {}

    # basis and potential
    s.basis = "DZVP-MOLOPT-SR-GTH"
    s.specific.cp2k.force_eval.dft.basis_set_file_name = "BASIS_MOLOPT"
    s.potential = "GTH-PBE"
    s.specific.cp2k.force_eval.dft.potential_file_name = "GTH_POTENTIALS"
    return s


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


def _read_result_file(result: Result, extension: str, max_line: int = 100) -> "None | str":
    """Find and read the first file in ``result`` with the provided file extension.

    Returns ``None`` if no such file can be found.
    """
    root = result.archive["plams_dir"]
    if root is None:
        return None

    iterator = (os.path.join(root, i) for i in os.listdir(root)
                if os.path.splitext(i)[1] == extension)
    for i in iterator:
        with open(i, "r") as f:
            ret_list = f.readlines()
            ret = "..." if len(ret_list) > max_line else ""
            ret += "".join(ret_list[-max_line:])
            return textwrap.indent(ret, 4 * " ")
    else:
        return None


def validate_status(result: Result, *, print_out: bool = True, print_err: bool = True) -> None:
    """Validate the status of the ``qmflows.Result`` object is set to ``"successful"``.

    Parameters
    ----------
    result : qmflows.Result
        The to-be validated ``Result`` object.
    print_out : bool
        Whether to included the content of the ``Result`` objects' .out file in the exception.
    print_err : bool
        Whether to included the content of the ``Result`` objects' .err file in the exception.

    Raises
    ------
    AssertionError
        Raised when :code:`result.status != "successful"`.

    """
    if result.status == "successful":
        return None

    indent = 4 * " "
    msg = f"Unexpected {result.job_name} status: {result.status!r}"

    if print_out:
        out = _read_result_file(result, ".out")
        if out is not None:
            msg += f"\n\nout_file:\n{out}"
    if print_err:
        err = _read_result_file(result, ".err")
        if err is not None:
            msg += f"\n\nerr_file:\n{err}"
    raise AssertionError(msg)


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
