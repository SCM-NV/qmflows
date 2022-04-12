"""Test package wrapper funcionality."""

from pathlib import Path

import pytest
from scm.plams import ADFJob, AMSJob, Molecule
from assertionlib import assertion

from qmflows import PackageWrapper, run, Settings
from qmflows.packages import ADF_Result, ResultWrapper
from qmflows.test_utils import PATH, requires_ams

WATER = Molecule(PATH / "water.xyz")
WATER.guess_bonds()


@pytest.mark.slow
@requires_ams
def test_package_wrapper(tmp_path: Path) -> None:
    """Tests for :class:`~qmflows.packages.PackageWrapper`."""
    s1 = Settings()
    s1.input.ams.Task = 'GeometryOptimization'
    s1.input.ams.GeometryOptimization.Convergence.Gradients = 1.0e-4
    s1.input.ams.Properties.NormalModes = 'true'
    s1.input.DFTB.Model = 'DFTB3'
    s1.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'

    s2 = Settings()
    s2.input.basis.type = 'DZP'
    s2.input.basis.core = 'None'
    s2.input.basis.createoutput = 'None'

    mol = WATER

    job1 = PackageWrapper(AMSJob)(s1, mol, name='amsjob')
    result1 = run(job1, path=tmp_path, folder="test_package_wrapper")
    assertion.isinstance(result1, ResultWrapper)
    assertion.eq(result1.results.job.status, 'successful')

    job2 = PackageWrapper(ADFJob)(s2, mol, name='adfjob')
    result2 = run(job2, path=tmp_path, folder="test_package_wrapper")
    assertion.isinstance(result2, ADF_Result)
    assertion.eq(result2.results.job.status, 'successful')
