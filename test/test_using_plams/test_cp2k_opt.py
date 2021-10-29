"""Run an small cp2k calculation."""

from pathlib import Path

from assertionlib import assertion
from scm.plams import Molecule

from qmflows import cp2k, run, templates
from qmflows.test_utils import PATH_MOLECULES, fill_cp2k_defaults, requires_cp2k


@requires_cp2k
def test_cp2k_opt(tmp_path: Path) -> None:
    """Run a simple molecular optimization."""
    s = fill_cp2k_defaults(templates.geometry)

    # Do a single step
    s.specific.cp2k.motion.geo_opt.max_iter = 1
    s.specific.cp2k.force_eval.dft.scf.eps_scf = 1e-1

    water = Molecule(PATH_MOLECULES / "h2o.xyz",
                     'xyz', charge=0, multiplicity=1)

    job = cp2k(s, water)
    mol = run(job, path=tmp_path, folder="test_cp2k_opt")
    assertion.isinstance(mol.molecule, Molecule)
