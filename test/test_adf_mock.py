"""Mock the ADF output."""
from assertionlib import assertion
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import adf, templates
from qmflows.packages.SCM import ADF_Result
from qmflows.test_utils import PATH, PATH_MOLECULES

WORKDIR = PATH / "output_adf"


def test_adf_mock(mocker):
    """Mock the ADF output."""
    mol = Molecule(PATH_MOLECULES / "acetonitrile.xyz")
    job = adf(templates.geometry, mol)

    run_mocked = mocker.patch("qmflows.run")
    jobname = "ADFjob"
    dill_path = WORKDIR / jobname / "ADFjob.dill"
    plams_dir = WORKDIR / jobname
    path_t21 = WORKDIR / jobname / "ADFjob.t21"
    run_mocked.return_value = ADF_Result(templates.geometry, mol, jobname, path_t21,
                                         dill_path=dill_path, plams_dir=plams_dir)
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

    # Test molecule
    # Molecule
    mol = rs.molecule
    assertion.isinstance(mol, Molecule)
    assertion.len_eq(mol, 6)  # there are 6 atoms
