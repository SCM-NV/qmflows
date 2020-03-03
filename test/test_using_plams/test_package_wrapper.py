import shutil
from os.path import isdir
from pathlib import Path

import pytest
from scm.plams import ADFJob, AMSJob, from_smiles

from qmflows import PackageWrapper, run, Settings
from qmflows.packages.SCM import ADF_Result
from qmflows.packages.package_wrapper import ResultWrapper

PATH = Path('test') / 'test_files'


@pytest.mark.slow
def test_package_wrapper() -> None:
    """Tests for :class:`PackageWrapper<qmflows.packages.package_wrapper.PackageWrapper>`."""
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

    mol = from_smiles('O')  # H2O

    workdir = PATH / 'workdir'

    try:
        job1 = PackageWrapper(AMSJob)(s1, mol, name='amsjob')
        result1 = run(job1, path=PATH, folder='workdir')

        assert isinstance(result1, ResultWrapper)
        assert result1.results.job.status == 'success'
    finally:
        if isdir(workdir):
            shutil.rmtree(workdir)

    try:
        job2 = PackageWrapper(ADFJob)(s2, mol, name='adfjob')
        result2 = run(job2, path=PATH, folder='workdir')

        assert isinstance(result2, ADF_Result)
        assert result2.results.job.status == 'success'
    finally:
        if isdir(workdir):
            shutil.rmtree(workdir)
