"""Check orca output reader."""
from assertionlib import assertion
from scm import plams

from qmflows import logger
from qmflows.parsers.generic_parsers import awk_file
from qmflows.parsers.orca_parser import (parse_hessian, parse_molecule,
                                         parse_molecule_traj)
from qmflows.test_utils import PATH

ORCA_WORKDIR = PATH / "output_orca" / "ORCAjob"
PATH_TRJ = ORCA_WORKDIR / "ORCAjob.trj"
PATH_OPT = ORCA_WORKDIR / "ORCAjob.opt"
PATH_OUTPUT = PATH / "output_orca" / "ORCAjob" / "ORCAjob.out"


def test_orca_mol_trj():
    """Test the result of reading a `job.trj` file."""
    mol = parse_molecule_traj(PATH_TRJ)
    logger.info(mol)

    assertion.isinstance(mol, plams.Molecule)


def test_orca_parse_hessian():
    """Check Hessian reader."""
    hess = parse_hessian(PATH_OPT, "$hessian_approx")
    assertion.shape_eq(hess, (18, 18))


def test_orca_awk():
    """Test the Awk reader for orca."""
    script = r"/Sum of individual times/ {runtime = $6} END {print runtime}"
    filename = PATH_OUTPUT.as_posix()
    result = awk_file(filename, script=script)
    assertion.isfinite(result)


def test_orca_parse_molecule():
    """Test the reading of the molecule from the output."""
    mol = parse_molecule(PATH_OUTPUT)
    assertion.isinstance(mol, plams.Molecule)
    assertion.len_eq(mol, 6)
