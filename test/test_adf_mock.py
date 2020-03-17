"""Mock the ADF output."""
from assertionlib import assertion
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import adf, templates
from qmflows.packages import Result, package_properties
from qmflows.test_utils import PATH, PATH_MOLECULES

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
    assertion.isfinite(rs.energy)
    assertion.isfinite(rs.homo)
    assertion.isfinite(rs.lumo)
    # Array of runtime for each optimization
    assertion.isfinite(rs.runtime[-1])
    # 3-dimensional vector
    assertion.len_eq(rs.dipole, 3)
    # number of steps until convergence
    assertion.eq(rs.optcycles, 8)
