"""Check orca output reader."""
from qmflows.parsers.orca_parser import (
    parse_hessian, parse_molecule_traj)
from scm import plams
from qmflows.test_utils import PATH

path_trj = PATH / "orca_output" / "ORCAjob.trj"
path_opt = PATH / "orca_output" / "ORCAjob.opt"


def test_orca_mol_trj():
    """Test the result of reading a `job.trj` file."""
    mol = parse_molecule_traj(path_trj)
    print(mol)

    assert isinstance(mol, plams.Molecule)


def test_parse_hessian():
    """Check Hessian reader."""
    hess = parse_hessian(path_opt, "$hessian_approx")
    assert hess.shape == (18, 18)
