"""Run an small cp2k calculation."""
from distutils.spawn import find_executable

import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import cp2k, run, templates
from qmflows.test_utils import PATH_MOLECULES, fill_cp2k_defaults, delete_output
from qmflows.type_hints import PathLike


def cp2k_available() -> None:
    """Check if cp2k is installed."""
    path = find_executable("cp2k.popt")

    return path is not None


@delete_output
@pytest.mark.skipif(
    not cp2k_available(), reason="CP2K is not install or not loaded")
def test_cp2k_opt(tmp_path: PathLike):
    """Run a simple molecular optimization."""
    s = fill_cp2k_defaults(templates.geometry)

    # Do a single step
    s.specific.cp2k.motion.geo_opt.max_iter = 1
    s.specific.cp2k.force_eval.dft.scf.eps_scf = 1e-1

    water = Molecule(PATH_MOLECULES / "h2o.xyz",
                     'xyz', charge=0, multiplicity=1)

    job = cp2k(s, water)
    mol = run(job, folder=tmp_path)
    assertion.isinstance(mol.molecule, Molecule)
