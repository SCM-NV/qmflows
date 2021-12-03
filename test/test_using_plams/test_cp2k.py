"""Run an small cp2k calculation."""

from pathlib import Path

import numpy as np
import h5py
import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, cp2k, run, templates
from qmflows.test_utils import PATH, PATH_MOLECULES, fill_cp2k_defaults, requires_cp2k


@requires_cp2k
@pytest.mark.slow
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


@requires_cp2k
@pytest.mark.slow
@pytest.mark.parametrize("mo_index_range", ["1_4", "1_2", "3_4"])
def test_cp2k_singlepoint(tmp_path: Path, mo_index_range: str) -> None:
    """Run a simple single point."""
    mol = Molecule(PATH_MOLECULES / "h2o.xyz", 'xyz', charge=0, multiplicity=1)

    s = fill_cp2k_defaults(templates.singlepoint)
    s.specific.cp2k.force_eval.dft.print.mo = Settings(
        add_last="numeric",
        each=Settings(qs_scf=0),
        eigenvalues="",
        eigenvectors="",
        filename="h2o",
        ndigits=36,
        occupation_numbers="",
        mo_index_range=mo_index_range.replace("_", " "),
    )
    job = cp2k(s, mol)
    result = run(job, path=tmp_path, folder="test_cp2k_singlepoint")

    with h5py.File(PATH / "test_output.hdf5", "r") as f:
        key = f"test_using_plams/test_cp2k/test_cp2k_singlepoint/{mo_index_range}"
        ref = f[key][...].view(np.recarray)

    orbitals = result.orbitals
    assertion.is_not(orbitals, None)
    np.testing.assert_allclose(orbitals.eigenvalues, ref.eigenvalues)
    np.testing.assert_allclose(np.abs(orbitals.eigenvectors).T, ref.eigenvectors)