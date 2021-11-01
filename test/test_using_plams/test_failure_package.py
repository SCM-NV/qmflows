"""Test that failure is handled correctly."""
import warnings
from pathlib import Path

import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, run, templates, logger
from qmflows.packages import adf, dftb, orca
from qmflows.test_utils import PATH_MOLECULES, requires_adf, requires_orca
from qmflows.warnings_qmflows import QMFlows_Warning


@pytest.mark.xfail
@requires_adf
def test_fail_scm(tmp_path: Path) -> None:
    """Test that both ADF and DFTB returns ``None`` if a computation fails."""
    # Temporary mute all QMFlows_Warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=QMFlows_Warning)

        mol = Molecule(PATH_MOLECULES / "ethylene.xyz")

        # Some dftb specific settings
        dftb_set = Settings()
        dftb_set.specific.dftb.dftb.scc

        # Calculate the DFTB hessian
        opt_dftb = dftb(templates.geometry.overlay(dftb_set), mol,
                        job_name="failed_DFTB")
        fail_adf = adf(None, opt_dftb.molecule, job_name="fail_adf")
        result = run(fail_adf.molecule, path=tmp_path, folder="test_fail_scm")
        logger.info(result)
        assertion.eq(result, None)


@requires_orca
def test_fail_orca(tmp_path: Path) -> None:
    """Orca package should returns ``None`` if the computation fails."""
    methanol = Molecule(PATH_MOLECULES / 'methanol.xyz')

    s = Settings()
    s.specific.orca.main = "RKS The_Cow_Functional SVP Opt TightSCF SmallPrint"

    opt = orca(s, methanol, job_name='fail_orca')

    result = run(opt.molecule, path=tmp_path, folder="test_fail_orca")
    assertion.eq(result, None)
