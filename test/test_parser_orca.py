from qmflows.parsers.orca_parser import (
    parse_hessian, parse_molecule_traj)
from scm import plams

path_trj = "test/test_files/ORCAjob.trj"
path_opt = "test/test_files/ORCAjob.opt"


def test_orca_mol_trj():
    """
    Test the result of reading a `job.trj` file.
    """
    mol = parse_molecule_traj(path_trj)
    print(mol)

    assert isinstance(mol, plams.Molecule)


def test_parse_hessian():

    hess = parse_hessian(path_opt, "$hessian_approx")
    assert hess.shape == (18, 18)
