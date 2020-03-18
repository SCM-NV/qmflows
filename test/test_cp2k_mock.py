"""Mock CP2K funcionality."""
import numpy as np
from assertionlib import assertion
from pytest_mock import mocker
from scm.plams import Molecule

from qmflows import cp2k, templates
from qmflows.packages import Result, package_properties
from qmflows.test_utils import PATH, PATH_MOLECULES
from qmflows.fileFunctions import yaml2Settings

# module constants
WORKDIR = PATH / "output_cp2k"

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


def test_cp2k_mock(mocker):
    """Mock a call to CP2K."""
    ETHYLENE = Molecule(PATH_MOLECULES / "ethylene.xyz")
    # geometry optimization input
    s = templates.singlepoint
    s.specific.cp2k.force_eval.subsys.cell.periodic = "None"
    s.cell_parameters = 10
    s = s.overlay(kinds_template)

    # print orbitals
    s.specific.cp2k.force_eval.dft.print.mo.filename = (
        PATH / "orbitals.out").as_posix()
    s.specific.cp2k.force_eval.dft.print.mo.mo_index_range = "7 46"
    s.specific.cp2k.force_eval.dft.scf.added_mos = 20

    # functional
    s.specific.cp2k.force_eval.dft.xc.xc_functional.pbe = {}

    # basis and potential
    s.specific.cp2k.force_eval.dft.potential_file_name = (
        PATH / "GTH_POTENTIALS").absolute().as_posix()
    s.specific.cp2k.force_eval.dft.basis_set_file_name = (
        PATH / "BASIS_MOLOPT").absolute().as_posix()

    job = cp2k(s, ETHYLENE)

    run_mocked = mocker.patch("qmflows.run")
    jobname = "cp2k_job"
    dill_path = WORKDIR / jobname / "cp2k_job.dill"
    plams_dir = WORKDIR / jobname
    adf_properties = package_properties["cp2k"]
    run_mocked.return_value = Result(templates.geometry, ETHYLENE, jobname, dill_path=dill_path,
                                     plams_dir=plams_dir, properties=adf_properties)
    rs = run_mocked(job)
    assertion.isfinite(rs.energy)

    # Molecular orbitals
    orbs = rs.orbitals
    assertion.isfinite(np.sum(orbs.eigenVals))  # eigenvalues
    assertion.shape_eq(orbs.coeffs, (46, 40))
