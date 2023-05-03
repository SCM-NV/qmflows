from pathlib import Path

import numpy as np
import pytest
from assertionlib import assertion
from scm.plams import Molecule

from qmflows import Settings, run, cp2k_mm, singlepoint, geometry, freq, md, cell_opt
from qmflows.test_utils import get_mm_settings, validate_status, PATH, PATH_MOLECULES, requires_cp2k

MOL = Molecule(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.xyz')

#: Example input Settings for CP2K mm calculations.
SETTINGS: Settings = get_mm_settings()


def overlap_coords(xyz1: np.ndarray, xyz2: np.ndarray) -> np.ndarray:
    """Rotate *xyz1* such that it overlaps with *xyz2* using the Kabsch algorithm."""
    xyz1 = np.asarray(xyz1, dtype=np.float64)
    xyz2 = np.asarray(xyz2, dtype=np.float64)

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
@requires_cp2k
def test_singlepoint(tmp_path: Path) -> None:
    """Test CP2K singlepoint calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += singlepoint.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL, job_name='cp2k_mm_sp')
    result = run(job, path=tmp_path, folder="test_singlepoint")
    validate_status(result)

    # Compare energies
    ref = -15.4431781758
    assertion.isclose(result.energy, ref, rel_tol=10**-4)


@pytest.mark.slow
@requires_cp2k
def test_geometry(tmp_path: Path) -> None:
    """Test CP2K geometry optimization calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += geometry.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL, job_name='cp2k_mm_opt')
    result = run(job, path=tmp_path, folder="test_geometry")
    validate_status(result)

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
@requires_cp2k
def test_freq(tmp_path: Path) -> None:
    """Test CP2K frequency calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += freq.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=mol, job_name='cp2k_mm_freq')
    result = run(job, path=tmp_path, folder="test_freq")
    validate_status(result)

    freqs = result.frequencies
    freqs_ref = np.load(PATH / 'Cd68Cl26Se55__26_acetate.freq.npy')
    np.testing.assert_allclose(freqs, freqs_ref, rtol=0, atol=5.0)

    G = result.free_energy
    H = result.enthalpy
    H_ref = -8642.371633064053
    assertion.isnan(G)
    assertion.isclose(H, H_ref, rel_tol=0.1)


@pytest.mark.slow
@requires_cp2k
def test_md(tmp_path: Path) -> None:
    """Test CP2K molecular dynamics calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(
        PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.freq.xyz')  # Optimized coordinates
    s = SETTINGS.copy()
    s.specific.cp2k += md.specific.cp2k_mm
    s.specific.cp2k.motion.md.steps = 1000

    job = cp2k_mm(settings=s, mol=mol, job_name='cp2k_mm_md')
    result = run(job, path=tmp_path, folder="test_md")
    validate_status(result)

    plams_results = result.results
    assertion.isfile(plams_results['cp2k-1_1000.restart'])


@pytest.mark.slow
@requires_cp2k
def test_c2pk_cell_opt(tmp_path: Path) -> None:
    """Test CP2K cell optimization calculations with the :class:`CP2K_MM` class."""
    mol = Molecule(PATH / 'cspbbr3_3d.xyz')

    s = Settings()
    s.specific.cp2k += cell_opt.specific.cp2k_mm.copy()
    s.specific.cp2k.motion.cell_opt.max_iter = 10
    s.specific.cp2k.motion.print['forces low'].filename = ''

    s.gmax = [22, 22, 22]
    s.cell_parameters = [25.452, 35.995, 24.452]
    s.charge = {
        'param': 'charge',
        'Cs': 0.2,
        'Pb': 0.4,
        'Br': -0.2,
    }
    s.lennard_jones = {
        'param': ('sigma', 'epsilon'),
        'unit': ('nm', 'kjmol'),
        'Cs Cs': (0.585, 1),
        'Cs Pb': (0.510, 1),
        'Br Se': (0.385, 1),
        'Pb Pb': (0.598, 1),
        'Br Pb': (0.290, 1),
        'Br Br': (0.426, 1),
    }

    job = cp2k_mm(settings=s, mol=mol, job_name='cp2k_mm_cell_opt')
    result = run(job, path=tmp_path, folder="test_c2pk_cell_opt")
    validate_status(result)
