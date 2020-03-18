"""Mock the DFTB output."""
import numpy as np
import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import dftb, templates
from qmflows.packages.SCM import DFTB_Result
from qmflows.test_utils import PATH

WORKDIR = PATH / "output_dftb"
WATER = molkit.from_smiles('[OH2]', forcefield='mmff')


def mock_runner(mocker_instance, jobname: str) -> DFTB_Result:
    """Create a Result instance using a mocked runner."""
    run_mocked = mocker_instance.patch("qmflows.run")
    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname
    run_mocked.return_value = DFTB_Result(templates.geometry, WATER, jobname, dill_path=dill_path,
                                          plams_dir=plams_dir)
    return run_mocked


def test_dftb_opt_mock(mocker):
    """Mock a geometry optimization using DFTB."""
    jobname = "dftb_geometry"
    job = dftb(templates.geometry, WATER, job_name=jobname)

    run_mocked = mock_runner(mocker, jobname)
    rs = run_mocked(job)

    # check the properties
    assertion.isfinite(rs.energy)
    assertion.len_eq(rs.dipole, 3)
    assertion.isinstance(rs.molecule, Molecule)


def test_dftb_freq_mock(mocker):
    """Mock a geometry optimization using DFTB."""
    jobname = "dftb_freq"
    job = dftb(templates.geometry, WATER, job_name=jobname)

    run_mocked = mock_runner(mocker, jobname)
    rs = run_mocked(job)
    # Check that hessian is symmetric
    hess_list = rs.hessian
    dim = int(np.sqrt(len(hess_list)))
    hess = np.array(hess_list).reshape(dim, dim)
    assertion.isclose(np.sum(hess - hess.T), 0.0)

    # enthalpy and Gibbs free energy
    assertion.isfinite(rs.enthalpy)
    assertion.isfinite(rs.free_energy)

    # frequencies
    assertion.isfinite(np.sum(rs.frequencies))
