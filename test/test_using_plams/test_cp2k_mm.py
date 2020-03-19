import numpy as np
import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, run, cp2k_mm, singlepoint, geometry, freq, md
from qmflows.test_utils import delete_output, get_mm_settings, PATH, PATH_MOLECULES

MOL = Molecule(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.xyz')

#: Example input Settings for CP2K mm calculations.
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


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_singlepoint() -> None:
    """Test CP2K singlepoint calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += singlepoint.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL, job_name='cp2k_mm_sp')
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Compare energies
    ref = -15.4431781758
    assertion.isclose(result.energy, ref, rel_tol=10**-4)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_geometry() -> None:
    """Test CP2K geometry optimization calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += geometry.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL, job_name='cp2k_mm_opt')
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Compare energies
    ref = -16.865587192150834
    assertion.isclose(result.energy, ref, rel_tol=10**-4)

    # Compare geometries
    xyz_ref = np.load(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.npy')
    _xyz = np.array(result.geometry)
    _xyz -= _xyz.mean(axis=0)[None, ...]
    xyz = overlap_coords(_xyz, xyz_ref)

    r_mean = np.linalg.norm(xyz - xyz_ref, axis=1).mean()
    assertion.le(r_mean, 0.1)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_freq() -> None:
    """Test CP2K frequency calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += freq.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=mol, job_name='cp2k_mm_freq')
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    freqs = result.frequencies
    freqs_ref = np.load(PATH / 'Cd68Cl26Se55__26_acetate.freq.npy')
    np.testing.assert_allclose(freqs, freqs_ref, rtol=0, atol=5.0)

    G = result.free_energy
    H = result.enthalpy
    H_ref = -8642.371633064053
    assertion.isnan(G)
    assertion.isclose(H, H_ref, rel_tol=0.1)


@pytest.mark.slow
@delete_output(delete_workdir=True)
def test_md() -> None:
    """Test CP2K frequency calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += md.specific.cp2k_mm
    s.specific.cp2k.motion.md.steps = 1000

    job = cp2k_mm(settings=s, mol=mol, job_name='cp2k_mm_md')
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    plams_results = result.results
    assertion.isfile(plams_results['cp2k-1_1000.restart'])
