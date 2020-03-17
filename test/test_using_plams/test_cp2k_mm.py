from os import PathLike
from io import TextIOBase
from itertools import islice
from typing import Union, AnyStr, Any
from collections import abc

import numpy as np
import pytest
from assertionlib import assertion
from scm.plams import Molecule, add_to_instance, add_to_class, Units, Cp2kJob

from qmflows import Settings, run, cp2k_mm, singlepoint, geometry, freq, md
from qmflows.backports import nullcontext
from qmflows.test_utils import delete_output, get_mm_settings, PATH, PATH_MOLECULES

MOL = Molecule(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.xyz')
SETTINGS: Settings = get_mm_settings()


def overlap_coords(xyz1: np.ndarray, xyz2: np.ndarray) -> np.ndarray:
    """Rotate *xyz1* such that it overlaps with *xyz2* using the Kabsch algorithm."""
    xyz1 = np.asarray(xyz1, dtype=float)
    xyz2 = np.asarray(xyz2, dtype=float)

    # Peform a singular value decomposition on the covariance matrix
    H = xyz1.T @ xyz2
    U, _, Vt = np.linalg.svd(H)
    V, Ut = Vt.T, U.T

    # Construct the rotation matrix
    rotmat = np.eye(3)
    rotmat[2, 2] = np.linalg.det(V @ Ut)
    rotmat = V @ rotmat @ Ut

    return xyz1 @ rotmat.T


@add_to_class(Cp2kJob)
def get_runscript(self):
    inp, out = self._filename('inp'), self._filename('out')
    return f'cp2k.ssmp -i {inp} -o {out}'


def get_energy(self, index: int = -1, unit: str = 'Hartree') -> float:
    """Return the energy of the last occurence of ``'ENERGY| Total FORCE_EVAL'`` in the output."""
    energy_str = self.grep_output('ENERGY| Total FORCE_EVAL')[index]
    energy = float(energy_str.rsplit(maxsplit=1)[1])
    return Units.convert(energy, 'Hartree', unit)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_singlepoint() -> None:
    """Test CP2K singlepoint calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += singlepoint.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Yes, this is a small hack as neither energy nor get_energy() seems to work
    plams_results = result.results
    add_to_instance(plams_results)(get_energy)
    setattr(result, 'get_energy', plams_results.get_energy)

    # Compare energies
    ref = -15.4431781758
    energy = result.get_energy()
    assertion.isclose(energy, ref, rel_tol=10**-4)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_geometry() -> None:
    """Test CP2K geometry optimization calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += geometry.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Yes, this is a small hack as neither energy nor get_energy() seems to work
    plams_results = result.results
    add_to_instance(plams_results)(get_energy)
    setattr(result, 'get_energy', plams_results.get_energy)

    # Compare energies
    ref = -16.865587192150834
    energy = result.get_energy()
    assertion.isclose(energy, ref, rel_tol=10**-4)

    # Compare geometries
    xyz_ref = np.load(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.npy')
    _xyz = Molecule(plams_results['cp2k-pos-1.xyz'], geometry=500).as_array()
    _xyz -= _xyz.mean(axis=0)[None, ...]
    xyz = overlap_coords(_xyz, xyz_ref)
    np.testing.assert_allclose(xyz, xyz_ref, rtol=np.inf, atol=0.05)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_freq() -> None:
    """Test CP2K frequency calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += freq.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=mol)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    plams_results = result.results

    freqs = get_frequencies(plams_results['cp2k-VIBRATIONS-1.mol'])
    freqs_ref = np.load(PATH / 'Cd68Cl26Se55__26_acetate.freq.npy')
    np.testing.assert_allclose(freqs, freqs_ref, rtol=np.inf, atol=5.0)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_md() -> None:
    """Test CP2K frequency calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += md.specific.cp2k_mm
    s.specific.cp2k.motion.md.steps = 1000

    job = cp2k_mm(settings=s, mol=mol)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    plams_results = result.results
    assertion.isfile(plams_results['cp2k-1_1000.restart'])
