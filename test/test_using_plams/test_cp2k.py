"""Run an small cp2k calculation."""

import warnings
from pathlib import Path

import pytest
import numpy as np
import h5py
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, cp2k, run, templates
from qmflows.warnings_qmflows import Orbital_Warning
from qmflows.test_utils import (
    PATH,
    PATH_MOLECULES,
    fill_cp2k_defaults,
    requires_cp2k,
    validate_status,
)


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
    result = run(job, path=tmp_path, folder="test_cp2k_opt")

    validate_status(result)
    assertion.isinstance(result.molecule, Molecule)


@requires_cp2k
@pytest.mark.slow
@pytest.mark.parametrize("mo_index_range", ["1_4", "1_2", "3_4", "0_99"])
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
    validate_status(result)

    with h5py.File(PATH / "test_output.hdf5", "r") as f:
        key = f"test_using_plams/test_cp2k/test_cp2k_singlepoint/{mo_index_range}"
        ref = f[key][...].view(np.recarray)

    if mo_index_range == "0_99":
        with pytest.warns(Orbital_Warning):
            orbitals = result.orbitals
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error", Orbital_Warning)
            orbitals = result.orbitals

    assertion.is_not(orbitals, None)
    np.testing.assert_allclose(orbitals.eigenvalues, ref.eigenvalues)
    np.testing.assert_allclose(np.abs(orbitals.eigenvectors).T, ref.eigenvectors)


@requires_cp2k
@pytest.mark.slow
def test_cp2k_freq(tmp_path: Path) -> None:
    """Run a simple single point."""
    mol = Molecule(PATH_MOLECULES / "h2o.xyz", 'xyz', charge=0, multiplicity=1)
    s = fill_cp2k_defaults(templates.freq)

    job = cp2k(s, mol)
    result = run(job, path=tmp_path, folder="test_cp2k_freq")

    validate_status(result)
    assertion.isclose(result.free_energy, -10801.971213467135)
    assertion.isclose(result.enthalpy, -10790.423489727531)
    np.testing.assert_allclose(result.frequencies, [1622.952012, 3366.668885, 3513.89377])
