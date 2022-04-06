"""Test Hessian utilities."""

from pathlib import Path

import pytest
import scm.plams.interfaces.molecule.rdkit as molkit
from noodles import gather
from scm.plams import Molecule

from qmflows import Settings, templates, logger
from qmflows.packages import run, orca, dftb
from qmflows.test_utils import PATH, requires_adf

WATER = Molecule(PATH / "water.xyz")
WATER.guess_bonds()


@pytest.mark.slow
@requires_adf
def test_hessian_transfer(tmp_path: Path) -> None:
    """Test DFTB -> Orca hessian transfer."""
    h2o = WATER.copy()
    h2o.properties.symmetry = 'C1'

    h2o_freq = dftb(templates.freq, h2o, job_name="freq").hessian

    s = Settings()
    s.inithess = h2o_freq

    h2o_opt = orca(templates.geometry.overlay(s), h2o, job_name="opt")

    energy = h2o_opt.energy
    dipole = h2o_opt.dipole

    wf = gather(energy, dipole)

    logger.info(run(wf, path=tmp_path, folder="test_hessian_transfer"))
