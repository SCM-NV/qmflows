"""Mock the ADF output."""
from scm.plams import Molecule

from pytest_mock import mocker
from qmflows import adf, templates
from qmflows.packages import Result, package_properties
from qmflows.test_utils import PATH, PATH_MOLECULES
import numpy as np

workdir = PATH / "output_adf"


def test_adf_mock(mocker):
    """Mock the ADF output."""
    mol = Molecule(PATH_MOLECULES / "acetonitrile.xyz")
    job = adf(templates.geometry, mol)

    run_mocked = mocker.patch("qmflows.run")
    jobname = "ADFjob"
    dill_path = workdir / jobname / "ADFjob.dill"
    plams_dir = workdir / jobname
    adf_properties = package_properties["adf"]
    run_mocked.return_value = Result(templates.geometry, mol, jobname, dill_path=dill_path,
                                     plams_dir=plams_dir, properties=adf_properties)
    rs = run_mocked(job)
    assert np.isfinite(rs.energy)
    assert np.isfinite(rs.homo)
    assert np.isfinite(rs.lumo)
    # Array of runtime for each optimization
    assert np.isfinite(rs.runtime[-1])
    # 3-dimensional vector
    assert len(rs.dipole) == 3
    # number of steps until convergence
    assert rs.optcycles == 8
