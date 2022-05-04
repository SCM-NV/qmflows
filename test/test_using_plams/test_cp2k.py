"""Run an small cp2k calculation."""

import sys
from pathlib import Path

import pytest
import numpy as np
import h5py
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, cp2k, run, templates
from qmflows.test_utils import (
    PATH,
    PATH_MOLECULES,
    fill_cp2k_defaults,
    requires_cp2k,
    validate_status,
)

if sys.version_info >= (3, 7):
    from builtins import dict as OrderedDict
else:
    from collections import OrderedDict


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

    orbitals = result.orbitals
    assertion.is_not(orbitals, None)
    np.testing.assert_allclose(orbitals.eigenvalues, ref.eigenvalues)
    np.testing.assert_allclose(np.abs(orbitals.eigenvectors).T, ref.eigenvectors)


UNRESTRICTED = OrderedDict(
    doublet=("doublet", 2, "ho"),
    triplet=("triplet", 3, "o2"),
)


@requires_cp2k
@pytest.mark.slow
@pytest.mark.parametrize("name,multiplicity,filename", UNRESTRICTED.values(), ids=UNRESTRICTED)
def test_cp2k_singlepoint_unrestricted(
    tmp_path: Path,
    name: str,
    multiplicity: int,
    filename: str,
) -> None:
    """Run a simple single point."""
    mol = Molecule(PATH_MOLECULES / f"{filename}.xyz", charge=0, multiplicity=multiplicity)

    s = fill_cp2k_defaults(templates.singlepoint)
    s.specific.cp2k.force_eval.dft.uks = ""
    s.specific.cp2k.force_eval.dft.multiplicity = multiplicity
    s.specific.cp2k.force_eval.dft.print.mo = Settings(
        add_last="numeric",
        each=Settings(qs_scf=0),
        eigenvalues="",
        eigenvectors="",
        filename=filename,
        ndigits=36,
        occupation_numbers="",
    )
    job = cp2k(s, mol)
    result = run(job, path=tmp_path, folder=f"test_cp2k_singlepoint_{name}")
    validate_status(result)

    orb_tup = result.orbitals
    assertion.is_not(orb_tup, None)

    with h5py.File(PATH / "test_output.hdf5", "r") as f:
        key = f"test_using_plams/test_cp2k/test_cp2k_singlepoint_unrestricted/{name}"
        ref_tup = (
            f[key + "/alpha"][...].view(np.recarray),
            f[key + "/beta"][...].view(np.recarray),
        )

    assertion.eq(len(orb_tup), len(ref_tup))
    for orbitals, ref in zip(orb_tup, ref_tup):
        np.testing.assert_allclose(orbitals.eigenvalues, ref.eigenvalues, rtol=0, atol=1e-03)
        np.testing.assert_allclose(
            np.abs(orbitals.eigenvectors).T, ref.eigenvectors, rtol=0, atol=1e-03
        )


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
