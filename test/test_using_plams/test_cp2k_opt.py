"""Run an small cp2k calculation."""
from distutils.spawn import find_executable

import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import cp2k, run, templates
from qmflows.test_utils import PATH, PATH_MOLECULES, fill_cp2k_defaults, delete_output
from qmflows.type_hints import PathLike
from qmflows.parsers.cp2KParser import get_cp2k_version_run


def cp2k_available() -> bool:
    """Check if cp2k is installed."""
    path = find_executable("cp2k.popt")
    return path is not None


HAS_CP2K = cp2k_available()
RUN_FILE = PATH / "output_cp2k" / "cp2k_freq" / "cp2k_freq.run"


@delete_output
@pytest.mark.skipif(not HAS_CP2K, reason="CP2K is not install or not loaded")
def test_cp2k_opt(tmp_path: PathLike) -> None:
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


@pytest.mark.skipif(not HAS_CP2K, reason="CP2K is not install or not loaded")
def test_get_cp2k_version() -> None:
    out = get_cp2k_version_run(RUN_FILE)
    assertion.truth(out)
