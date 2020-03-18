"""Mock the DFTB output."""
import scm.plams.interfaces.molecule.rdkit as molkit
from assertionlib import assertion
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import dftb, templates
from qmflows.packages.SCM import DFTB_Result
from qmflows.test_utils import PATH, delete_output

WORKDIR = PATH / "output_dftb"


@delete_output
def test_dftb_opt_mock(mocker):
    """Mock a geometry optimization using DFTB."""
    mol = molkit.from_smiles('[OH2]', forcefield='mmff')
    jobname = "dftb_geometry"
    job = dftb(templates.geometry, mol, job_name=jobname)

    run_mocked = mocker.patch("qmflows.run")
    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname
    run_mocked.return_value = DFTB_Result(templates.geometry, mol, jobname,
                                          dill_path=dill_path, plams_dir=plams_dir)
    rs = run_mocked(job)
    assertion.isfinite(rs.energy)
    assertion.len_eq(rs.dipole, 3)
    assertion.isinstance(rs.molecule, Molecule)