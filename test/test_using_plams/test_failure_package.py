"""Test that failure is handled correctly."""
import warnings

import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, run, templates
from qmflows.packages import adf, dftb, orca
from qmflows.test_utils import PATH_MOLECULES, delete_output
from qmflows.warnings_qmflows import QMFlows_Warning


@delete_output
@pytest.mark.xfail
def test_fail_scm(tmpdir):
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
        result = run(fail_adf.molecule, path=tmpdir)
        print(result)
        assertion.eq(result, None)


@delete_output
def test_fail_orca(tmpdir):
    """Orca package should returns ``None`` if the computation fails."""
    methanol = Molecule(PATH_MOLECULES / 'methanol.xyz')

    s = Settings()
    s.specific.orca.main = "RKS The_Cow_Functional SVP Opt TightSCF SmallPrint"

    opt = orca(s, methanol, job_name='fail_orca')

    result = run(opt.molecule, path=tmpdir)
    assertion.eq(result, None)
