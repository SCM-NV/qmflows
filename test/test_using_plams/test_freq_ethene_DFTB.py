"""Test freq calculations."""

import pytest
from scm.plams import Molecule
from assertionlib import assertion
from qmflows import (dftb, run, templates)
from qmflows.test_utils import PATH_MOLECULES


@pytest.mark.slow
def test_freq():
    """Do some constraint optimizations then launch a freq calc."""
    mol = Molecule(PATH_MOLECULES / "ethene.xyz", "xyz")
    geo_opt = dftb(templates.geometry, mol)
    freq_calc = dftb(templates.freq, geo_opt.molecule, job_name="freq")
    r = run(freq_calc)
    assertion.eq(int(r.frequencies[0]), 831)
    assertion.len_eq(r.frequencies, 12)
