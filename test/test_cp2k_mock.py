"""Mock CP2K funcionality."""
import copy
from typing import Callable

import pytest
import numpy as np
from assertionlib import assertion
from pytest_mock import MockFixture
from scm.plams import Molecule

from qmflows import cp2k, templates
from qmflows.packages.cp2k_package import CP2K_Result
from qmflows.utils import init_restart
from qmflows.common import CP2KVersion
from qmflows.parsers.cp2KParser import get_cp2k_version_run
from qmflows.test_utils import (
    PATH,
    PATH_MOLECULES,
    fill_cp2k_defaults,
    requires_cp2k,
    validate_status,
)

try:
    CP2K_VERSION = get_cp2k_version_run("cp2k.popt")
except Exception:
    CP2K_VERSION = CP2KVersion(0, 0)


def mock_runner(mocker_instance, jobname: str) -> Callable[..., CP2K_Result]:
    """Create a Result instance using a mocked runner."""
    run_mocked = mocker_instance.patch("qmflows.run")
    run_mocked.return_value = CP2K_Result(
        templates.geometry,
        ETHYLENE,
        jobname,
        status="successful",
        dill_path=WORKDIR / jobname / f"{jobname}.dill",
        plams_dir=WORKDIR / jobname,
    )
    return run_mocked

# module constants
WORKDIR = PATH / "output_cp2k"
ETHYLENE = Molecule(PATH_MOLECULES / "ethylene.xyz")


def test_deepcopy():
    """Test the copy of a result instance."""
    jobname = "cp2k_job"
    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname
    result = CP2K_Result(templates.geometry, ETHYLENE, jobname,
                         dill_path=dill_path, plams_dir=plams_dir)

    copy_result = copy.deepcopy(result)

    init_restart(folder=WORKDIR)
    assertion.is_(copy_result.molecule, result.molecule)
    assertion.eq(copy_result.archive, result.archive)
    assertion.eq(copy_result.settings, result.settings)


# See https://github.com/SCM-NV/qmflows/issues/225
@pytest.mark.xfail(raises=AssertionError)
@pytest.mark.skipif(CP2K_VERSION >= (8, 2), reason="Requires CP2K < 8.2")
@requires_cp2k
def test_cp2k_singlepoint_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # single point calculation
    s = fill_cp2k_defaults(templates.singlepoint)

    # print orbitals
    s.specific.cp2k.force_eval.dft.print.mo.filename = (
        PATH / "orbitals.out").as_posix()
    s.specific.cp2k.force_eval.dft.print.mo.mo_index_range = "7 46"
    s.specific.cp2k.force_eval.dft.scf.added_mos = 20

    # Construct a result objects
    job = cp2k(s, ETHYLENE)
    jobname = "cp2k_job"
    run_mocked = mock_runner(mocker, jobname)
    rs = run_mocked(job)
    validate_status(rs)

    # electronic energy
    assertion.isfinite(rs.energy)

    # Molecular orbitals
    orbs = rs.orbitals
    assertion.assert_(np.isfinite, orbs.eigenvectors, post_process=np.all)  # eigenvalues
    assertion.shape_eq(orbs.eigenvectors, (46, 40))


def test_c2pk_opt_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # geometry optimization input
    s = fill_cp2k_defaults(templates.geometry)

    jobname = "cp2k_opt"
    run_mocked = mock_runner(mocker, jobname)

    job = cp2k(s, ETHYLENE, job_name=jobname)
    rs = run_mocked(job)
    validate_status(rs)

    # check the optimized geometry
    mol = rs.geometry
    assertion.len_eq(mol, 6)
    atom = mol[1]
    assertion.len_eq(atom.coords, 3)
    assertion.eq(atom.symbol, 'C')


def test_c2pk_freq_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # Frequency calculation
    s = fill_cp2k_defaults(templates.singlepoint)
    s.specific.cp2k.vibrational_analysis.thermochemistry = ".TRUE."
    s.specific.cp2k["global"]["run_type"] = "VIBRATIONAL_ANALYSIS"

    jobname = "cp2k_freq"
    run_mocked = mock_runner(mocker, jobname)
    job = cp2k(s, ETHYLENE, job_name=jobname)
    rs = run_mocked(job)
    validate_status(rs)

    # check properties
    assertion.isfinite(rs.enthalpy)
    assertion.isfinite(rs.free_energy)


def test_dir(mocker: MockFixture) -> None:
    s = fill_cp2k_defaults(templates.geometry)
    jobname = "cp2k_opt"
    run_mocked = mock_runner(mocker, jobname)

    job = cp2k(s, ETHYLENE, job_name=jobname)
    r = run_mocked(job)
    validate_status(r)

    assertion.issubset(CP2K_Result.prop_mapping.keys(), dir(r))
