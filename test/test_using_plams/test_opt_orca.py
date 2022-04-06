"""Tests for Orca functionality."""

from math import sqrt
from pathlib import Path

import pytest
from more_itertools import collapse
from scm.plams import Molecule
from assertionlib import assertion

from qmflows import Settings, templates, logger
from qmflows.packages import run, orca, dftb
from qmflows.test_utils import PATH_MOLECULES, requires_orca


@pytest.mark.slow
@requires_orca
def test_opt_orca(tmp_path: Path) -> None:
    """Test Orca input generation and run functions."""
    h2o = Molecule(PATH_MOLECULES / "h2o.xyz",
                   'xyz', charge=0, multiplicity=1)

    h2o_geometry = dftb(templates.geometry, h2o)

    s = Settings()
    # generic keyword "basis" must be present in the generic dictionary
    s.basis = "sto_dzp"
    # s.specific.adf.basis.core = "large"

    r = templates.singlepoint.overlay(s)
    h2o_singlepoint = orca(r, h2o_geometry.molecule)

    dipole = h2o_singlepoint.dipole

    final_result = run(dipole, n_processes=1, path=tmp_path, folder="test_opt_orca")

    expected_dipole = [0.82409, 0.1933, -0.08316]
    diff = sqrt(sum((x - y) ** 2 for x, y in zip(final_result,
                                                 expected_dipole)))
    logger.info(f"Expected dipole computed with Orca 3.0.3 is: {expected_dipole}")
    logger.info(f"Actual dipole is: {final_result}")

    assertion.lt(diff, 1e-2)


@pytest.mark.slow
@requires_orca
def test_methanol_opt_orca(tmp_path: Path) -> None:
    """Run a methanol optimization and retrieve the optimized geom."""
    methanol = Molecule(PATH_MOLECULES / "methanol.xyz")

    s = Settings()
    s.specific.orca.main = "RKS B3LYP SVP Opt NumFreq TightSCF SmallPrint"
    s.specific.orca.scf = " print[p_mos] 1"

    opt = orca(s, methanol)

    # extract coordinates
    mol_opt = run(opt.molecule, path=tmp_path, folder="test_methanol_opt_orca")
    coords = collapse([a.coords for a in mol_opt.atoms])
    logger.info(coords)
