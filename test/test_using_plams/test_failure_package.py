from scm.plams import Molecule
from qmflows import (Settings, templates, run)
from qmflows.packages import (adf, dirac, dftb, gamess, orca)
import os
import pytest
import shutil
import tempfile


def isNone(x):
    return True if x is None else False


def remove(folder):
    if os.path.isdir(folder):
        shutil.rmtree(folder)


@pytest.mark.xfail
def test_fail_scm():
    """ Test that both ADF and DFTB returns ``None`` if a computation fails"""
    # 5 membered ring from which ozone will dissociate
    folder = tempfile.mkdtemp(prefix="qmflows_")
    mol = Molecule("test/test_files/ethylene.xyz")

    # Some dftb specific settings
    dftb_set = Settings()
    dftb_set.specific.dftb.dftb.scc

    # Calculate the DFTB hessian
    opt_dftb = dftb(templates.geometry.overlay(dftb_set), mol,
                    job_name="failed_DFTB")
    fail_adf = adf(None, opt_dftb.molecule, job_name="fail_adf")
    try:
        result = run(fail_adf.molecule, path=folder)
        print(result)
        assert isNone(result)
    finally:
        remove(folder)


@pytest.mark.xfail
def test_fail_dirac():
    """ Dirac package should return ``None`` if it fails """
    folder = tempfile.mkdtemp(prefix="qmflows_")
    mol = Molecule("test/test_files/h2.xyz")
    s = Settings()
    job = dirac(s, mol, job_name="fail_dirac")
    try:
        result = run(job.energy, path=folder)
        assert isNone(result)
    finally:
        remove(folder)


def test_fail_gamess():
    """
    Gamess should return ``None`` if a calculation fails.
    """
    folder = tempfile.mkdtemp(prefix="qmflows_")
    symmetry = "Cpi"  # Erroneous Keyowrkd
    methanol = Molecule('test/test_files/ion_methanol.xyz')
    methanol.properties['symmetry'] = symmetry

    s = Settings()
    s.specific.gamess.contrl.nzvar = 12
    s.specific.gamess.pcm.solvnt = 'water'
    s.specific.gamess.basis.gbasis = 'sto'
    s.specific.gamess.basis.ngauss = 3

    inp = templates.geometry.overlay(s)
    methanol_geometry = gamess(inp, methanol, job_name="fail_gamess",
                               work_dir='/tmp')

    try:
        result = run(methanol_geometry.molecule, path=folder)
        assert isNone(result)
    finally:
        remove(folder)


def test_fail_orca():
    """ Orca package should returns ``None`` if the computation fails"""
    folder = tempfile.mkdtemp(prefix="qmflows_")

    methanol = Molecule('test/test_files/methanol.xyz')

    s = Settings()
    s.specific.orca.main = "RKS The_Cow_Functional SVP Opt TightSCF SmallPrint"

    opt = orca(s, methanol, job_name='fail_orca')

    try:
        result = run(opt.molecule, path=folder)
        assert isNone(result)
    finally:
        remove(folder)


if __name__ == "__main__":
    test_fail_gamess()
