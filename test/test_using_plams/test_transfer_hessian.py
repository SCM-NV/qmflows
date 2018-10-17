from qmflows import (Settings, templates)
from noodles import gather

# User Defined imports
from qmflows.packages.SCM import dftb
from qmflows.packages.orca import orca
from qmflows.packages import run
import pytest
import scm.plams.interfaces.molecule.rdkit as molkit


@pytest.mark.slow
def test_hessian_transfer():
    """
    Test DFTB -> Orca hessian transfer
    """
    h2o = molkit.from_smiles('O')
    h2o.properties.symmetry = 'C1'

    h2o_freq = dftb(templates.freq, h2o, job_name="freq").hessian

    s = Settings()
    s.inithess = h2o_freq

    h2o_opt = orca(templates.geometry.overlay(s), h2o, job_name="opt")

    energy = h2o_opt.energy
    dipole = h2o_opt.dipole

    wf = gather(energy, dipole)

    print(run(wf))
