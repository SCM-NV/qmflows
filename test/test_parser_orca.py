"""Check orca output reader."""
import numpy as np
from scm import plams

from qmflows.parsers.generic_parsers import awk_file
from qmflows.parsers.orca_parser import parse_hessian, parse_molecule_traj
from qmflows.test_utils import PATH

ORCA_WORKDIR = PATH / "output_orca" / "ORCAjob"
path_trj = ORCA_WORKDIR / "ORCAjob.trj"
path_opt = ORCA_WORKDIR / "ORCAjob.opt"


def test_orca_mol_trj():
    """Test the result of reading a `job.trj` file."""
    mol = parse_molecule_traj(path_trj)
    print(mol)

    assert isinstance(mol, plams.Molecule)


def test_parse_hessian():
    """Check Hessian reader."""
    hess = parse_hessian(path_opt, "$hessian_approx")
    assert hess.shape == (18, 18)


def test_awk_orca():
    """Test the Awk reader for orca."""
    script = r"/Sum of individual times/ {runtime = $6} END {print runtime}"
    path_output = PATH / "output_orca" / "ORCAjob"
    filename = (path_output / "ORCAjob.out").as_posix()
    result = awk_file(filename, script=script)
    assert np.isfinite(result)
