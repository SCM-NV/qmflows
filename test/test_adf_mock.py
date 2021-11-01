"""Mock the ADF output."""

from pathlib import Path
from typing import Optional

import pytest
from assertionlib import assertion
from pytest_mock import MockFixture
from scm.plams import Molecule

from qmflows import adf, templates
from qmflows.packages.SCM import ADF_Result
from qmflows.test_utils import PATH, PATH_MOLECULES
from qmflows.warnings_qmflows import QMFlows_Warning
from qmflows.utils import InitRestart

WORKDIR = PATH / "output_adf"


def test_adf_mock(mocker: MockFixture):
    """Mock the ADF output."""
    mol = Molecule(PATH_MOLECULES / "acetonitrile.xyz")
    job = adf(templates.geometry, mol)

    run_mocked = mocker.patch("qmflows.run")
    jobname = "ADFjob"
    dill_path = WORKDIR / jobname / "ADFjob.dill"
    plams_dir = WORKDIR / jobname
    run_mocked.return_value = ADF_Result(templates.geometry, mol, jobname,
                                         dill_path=dill_path, work_dir=plams_dir,
                                         plams_dir=plams_dir)
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


@pytest.mark.parametrize("name,match,status", [
    ("bob", "Generic property 'bob' not defined", None),
    ("energy", "It is not possible to retrieve property: 'energy'", "crashed"),
], ids=["undefined_property", "job_crashed"])
def test_getattr_warning(tmp_path: Path, name: str, match: str, status: Optional[str]) -> None:
    mol = Molecule(PATH_MOLECULES / "acetonitrile.xyz")
    jobname = "ADFjob"
    dill_path = WORKDIR / jobname / "ADFjob.dill"
    plams_dir = WORKDIR / jobname
    result = ADF_Result(templates.geometry, mol, jobname,
                        dill_path=dill_path, work_dir=plams_dir,
                        plams_dir=plams_dir)
    if status is not None:
        result.status = status

    # Need to fire up `plams.init` so the .dill file can be unpickled by `plams.load_job`
    with InitRestart(path=tmp_path, folder="test_getattr_warning"):
        with pytest.warns(QMFlows_Warning, match=match) as rec:
            getattr(result, name)
            assertion.len_eq(rec.list, 1)
            assertion.eq(rec.list[0].filename, __file__)

        # Test a second time after the .dill file has been unpickled
        with pytest.warns(QMFlows_Warning, match=match) as rec:
            getattr(result, name)
            assertion.len_eq(rec.list, 1)
            assertion.eq(rec.list[0].filename, __file__)
