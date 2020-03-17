"""Mock orca funcionality."""
import numpy as np
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import Settings, orca
from qmflows.packages import Result, package_properties
from qmflows.test_utils import PATH, PATH_MOLECULES

workdir = PATH / "output_orca"


def test_orca_mock(mocker):
    """Mock a call to orca."""
    methanol = Molecule(PATH_MOLECULES / "methanol.xyz")

    s = Settings()
    s.specific.orca.main = "RKS B3LYP SVP Opt NumFreq TightSCF SmallPrint"
    # print the orbitals
    s.specific.orca.scf = " print[p_mos] 1"
    job = orca(s, methanol)

    run_mocked = mocker.patch("qmflows.run")
    jobname = "ORCAjob"
    dill_path = workdir / jobname / "ORCAjob.dill"
    plams_dir = workdir / jobname
    adf_properties = package_properties["orca"]
    run_mocked.return_value = Result(s, methanol, jobname, dill_path=dill_path,
                                     plams_dir=plams_dir, properties=adf_properties)
    rs = run_mocked(job)

    assert np.isfinite(rs.energy)
    assert np.isfinite(rs.runtime)
    assert np.isfinite(np.sum(rs.dipole))
    # steps until convergence
    assert rs.optcycles == 8

    # Check that hessian is symmetric
    hess = rs.hessian
    assert np.isclose(np.sum(hess - hess.T), 0.0)

    # frequencies
    frequencies = rs.frequencies
    assert len(frequencies) == 18

    # Normal modes
    normal_modes = rs.normal_modes
    assert normal_modes.shape == (18, 18)

    # Orbitals
    orbs = rs.orbitals
    assert np.isfinite(np.sum(orbs.eigenVals))  # eigenvalues
