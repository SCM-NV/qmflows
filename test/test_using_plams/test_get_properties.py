"""Test the rdkit/Molecule interface."""

from pathlib import Path

import pytest
from scm.plams import Molecule

from qmflows import adf, dftb, run, templates
from qmflows.test_utils import PATH, requires_adf

HF = Molecule(PATH / "HF.xyz")
HF.guess_bonds()


@pytest.mark.slow
@requires_adf
def test_dftb_props(tmp_path: Path) -> None:
    """Get properties from DFTB freq calc."""
    mol = HF.copy()
    job = dftb(templates.freq, mol, job_name='dftb_FH')
    result = run(job, path=tmp_path, folder="test_dftb_props")

    expected_energy = -4.76
    assert abs(result.energy - expected_energy) < 0.01
    assert len(result.dipole) == 3
    expected_frequency = 3460.92
    assert abs(result.frequencies - expected_frequency) < 0.1
    assert len(result.charges) == 2
    assert abs(result.charges[0] + result.charges[1]) < 1e-6


@pytest.mark.slow
@requires_adf
def test_adf_props(tmp_path: Path) -> None:
    """Get properties from ADF freq calc."""
    mol = HF.copy()
    job = adf(templates.freq, mol, job_name='adf_FH')
    result = run(job, path=tmp_path, folder="test_adf_props")
    expected_energy = -0.30
    assert abs(result.energy - expected_energy) < 0.01
    assert len(result.dipole) == 3
    expected_frequency = 3480.90
    assert abs(result.frequencies[1] - expected_frequency) < 0.1
    assert len(result.charges) == 2
    assert abs(result.charges[0] + result.charges[1]) < 1e-6
