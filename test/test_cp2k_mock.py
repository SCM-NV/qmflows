"""Mock CP2K funcionality."""
import copy

import numpy as np
from assertionlib import assertion
from pytest_mock import MockFixture, mocker
from scm.plams import Molecule

from qmflows import Settings, cp2k, templates
from qmflows.fileFunctions import yaml2Settings
from qmflows.packages.cp2k_package import CP2K_Result
from qmflows.test_utils import PATH, PATH_MOLECULES
from qmflows.utils import init_restart

# module constants
WORKDIR = PATH / "output_cp2k"
ETHYLENE = Molecule(PATH_MOLECULES / "ethylene.xyz")

kinds_template = yaml2Settings("""
specific:
  cp2k:
     force_eval:
       subsys:
         kind:
           C:
             basis_set: DZVP-MOLOPT-SR-GTH-q4
             potential: GTH-PBE-q4
           H:
             basis_set: DZVP-MOLOPT-SR-GTH-q1
             potential: GTH-PBE-q1
""")


def fill_defaults(s: Settings) -> Settings:
    """Fill missing values from a job template."""
    s.periodic = "None"
    s.cell_parameters = 10
    s = s.overlay(kinds_template)

    # functional
    s.specific.cp2k.force_eval.dft.xc.xc_functional.pbe = {}

    # basis and potential
    add_basis_potential(s)

    return s


def add_basis_potential(s: Settings) -> None:
    """Add basis and potential path to the settings."""
    s.specific.cp2k.force_eval.dft.potential_file_name = (
        PATH / "GTH_POTENTIALS").absolute().as_posix()
    s.specific.cp2k.force_eval.dft.basis_set_file_name = (
        PATH / "BASIS_MOLOPT").absolute().as_posix()


def mock_runner(mocker_instance, jobname: str) -> CP2K_Result:
    """Create a Result instance using a mocked runner."""
    run_mocked = mocker_instance.patch("qmflows.run")
    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname
    run_mocked.return_value = CP2K_Result(templates.geometry, ETHYLENE, jobname,
                                          dill_path=dill_path, plams_dir=plams_dir)

    return run_mocked


def test_deepcopy():
    """Test the copy of a result instance."""
    jobname = "cp2k_job"
    dill_path = WORKDIR / jobname / f"{jobname}.dill"
    plams_dir = WORKDIR / jobname
    result = CP2K_Result(templates.geometry, ETHYLENE, jobname,
                         dill_path=dill_path, plams_dir=plams_dir)

    copy_result = copy.deepcopy(result)

    init_restart(folder=WORKDIR)
    assertion.is_(copy_result.molecule, result.molecule)
    assertion.eq(copy_result.archive, result.archive)
    assertion.eq(copy_result.settings, result.settings)


def test_cp2k_singlepoint_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # single point calculation
    s = fill_defaults(templates.singlepoint)

    # print orbitals
    s.specific.cp2k.force_eval.dft.print.mo.filename = (
        PATH / "orbitals.out").as_posix()
    s.specific.cp2k.force_eval.dft.print.mo.mo_index_range = "7 46"
    s.specific.cp2k.force_eval.dft.scf.added_mos = 20

    job = cp2k(s, ETHYLENE)
    jobname = "cp2k_job"
    run_mocked = mock_runner(mocker, jobname)
    rs = run_mocked(job)

    # electronic energy
    assertion.isfinite(rs.energy)

    # Molecular orbitals
    orbs = rs.orbitals
    assertion.isfinite(np.sum(orbs.eigenVals))  # eigenvalues
    assertion.shape_eq(orbs.coeffs, (46, 40))


def test_c2pk_opt_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # geometry optimization input
    s = fill_defaults(templates.geometry)

    jobname = "cp2k_opt"
    run_mocked = mock_runner(mocker, jobname)

    job = cp2k(s, ETHYLENE, job_name=jobname)
    rs = run_mocked(job)
    # check the optimized geometry
    mol = rs.geometry
    assertion.len_eq(mol, 6)
    atom = mol[1]
    assertion.len_eq(atom.coords, 3)
    assertion.eq(atom.symbol, 'C')


def test_c2pk_freq_mock(mocker: MockFixture):
    """Mock a call to CP2K."""
    # Frequency calculation
    s = fill_defaults(templates.singlepoint)
    s.specific.cp2k.vibrational_analysis.thermochemistry = ".TRUE."
    s.specific.cp2k["global"]["run_type"] = "VIBRATIONAL_ANALYSIS"

    jobname = "cp2k_freq"
    run_mocked = mock_runner(mocker, jobname)
    job = cp2k(s, ETHYLENE, job_name=jobname)
    rs = run_mocked(job)

    # check properties
    assertion.isfinite(rs.enthalpy)
    assertion.isfinite(rs.free_energy)
