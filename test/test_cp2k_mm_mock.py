"""Mock CP2K funcionality."""

import os
import shutil
from typing import Callable

import numpy as np
from assertionlib import assertion
from pytest_mock import MockFixture
from scm.plams import Molecule

from qmflows import Settings, cp2k_mm, singlepoint, geometry, freq, md, cell_opt
from qmflows.utils import InitRestart
from qmflows.packages.cp2k_mm import CP2KMM_Result
from qmflows.test_utils import get_mm_settings, validate_status, PATH, PATH_MOLECULES

MOL = Molecule(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.xyz')
WORKDIR = PATH / 'output_cp2k_mm'

#: Example input Settings for CP2K mm calculations.
SETTINGS: Settings = get_mm_settings()

# Ensure that plams.config is populated with a JobManager
with InitRestart(PATH, 'tmp'):
    pass
if os.path.isdir(PATH / 'tmp'):
    shutil.rmtree(PATH / 'tmp')


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


def mock_runner(mocker_instance: MockFixture,
                settings: Settings = SETTINGS,
                jobname: str = 'job') -> Callable[..., CP2KMM_Result]:
    """Create a Result instance using a mocked runner."""
    run_mocked = mocker_instance.patch("qmflows.run")

    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname

    run_mocked.return_value = CP2KMM_Result(
        settings, MOL, jobname,
        dill_path=dill_path,
        plams_dir=plams_dir,
        work_dir=WORKDIR,
        status='successful'
    )
    return run_mocked


def test_cp2k_singlepoint_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
    s = SETTINGS.copy()
    s.specific.cp2k += geometry.specific.cp2k_mm

    job = cp2k_mm(s, MOL)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_sp")

    result = run_mocked(job)
    validate_status(result)

    # Compare energies
    ref = -15.4431781758
    assertion.isclose(result.energy, ref, rel_tol=10**-4)


def test_c2pk_opt_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
    s = SETTINGS.copy()
    s.specific.cp2k += singlepoint.specific.cp2k_mm

    job = cp2k_mm(s, MOL)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_opt")

    result = run_mocked(job)
    validate_status(result)

    # Compare energies
    ref = -16.865587192150834
    energy = result.energy
    assertion.isclose(energy, ref, rel_tol=10**-4)

    # Compare geometries
    xyz_ref = np.load(PATH_MOLECULES / 'Cd68Cl26Se55__26_acetate.npy')
    _xyz = np.array(result.geometry)
    _xyz -= _xyz.mean(axis=0)[None, ...]
    xyz = overlap_coords(_xyz, xyz_ref)

    r_mean = np.linalg.norm(xyz - xyz_ref, axis=1).mean()
    assertion.le(r_mean, 0.1)


def test_c2pk_freq_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
    s = SETTINGS.copy()
    s.specific.cp2k += freq.specific.cp2k_mm

    job = cp2k_mm(s, MOL)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_freq")

    result = run_mocked(job)
    validate_status(result)

    freqs = result.frequencies
    freqs_ref = np.load(PATH / 'Cd68Cl26Se55__26_acetate.freq.npy')
    np.testing.assert_allclose(freqs, freqs_ref, rtol=0, atol=5.0)

    G = result.free_energy
    H = result.enthalpy
    H_ref = -8642.371633064053
    assertion.isnan(G)
    assertion.isclose(H, H_ref, rel_tol=0.1)


def test_c2pk_md_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
    s = SETTINGS.copy()
    s.specific.cp2k += md.specific.cp2k_mm

    job = cp2k_mm(s, MOL)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_md")

    result = run_mocked(job)
    validate_status(result)

    assertion.isfile(result.results['cp2k-1_1000.restart'])


def test_c2pk_cell_opt_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
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

    job = cp2k_mm(s, mol)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_cell_opt")
    result = run_mocked(job)
    validate_status(result)

    ref_volume = np.load(PATH / 'volume.npy')
    ref_coordinates = np.load(PATH / 'coordinates.npy')
    ref_forces = np.load(PATH / 'forces.npy')
    ref_lattice = np.load(PATH / 'lattice.npy')

    np.testing.assert_allclose(result.volume, ref_volume)
    np.testing.assert_allclose(result.coordinates, ref_coordinates)
    np.testing.assert_allclose(result.forces, ref_forces)
    np.testing.assert_allclose(result.lattice, ref_lattice)


def test_c2pk_npt_mock(mocker: MockFixture) -> None:
    """Mock a call to CP2K."""
    s = Settings()
    job = cp2k_mm(s, None)
    run_mocked = mock_runner(mocker, settings=s, jobname="cp2k_mm_npt")
    result = run_mocked(job)
    validate_status(result)

    ref_pressure = np.load(PATH / 'pressure.npy')
    np.testing.assert_allclose(result.pressure, ref_pressure)
