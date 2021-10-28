"""Test Hessian utilities."""
import pytest
import scm.plams.interfaces.molecule.rdkit as molkit
from noodles import gather
from scm.plams import Molecule

from qmflows import Settings, templates, logger
from qmflows.packages import run
from qmflows.packages.orca import orca
from qmflows.packages.SCM import dftb
from qmflows.test_utils import PATH, requires_adf

WATER = Molecule(PATH / "water.xyz")
WATER.guess_bonds()


@pytest.mark.slow
@requires_adf
def test_hessian_transfer():
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

    logger.info(run(wf))
