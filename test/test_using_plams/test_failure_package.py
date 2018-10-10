from scm.plams import Molecule
from qmflows import (Settings, templates, run)

# User Defined imports
from qmflows.packages import (adf, dirac, dftb, gamess, orca)


def isNone(x):
    return True if x is None else False


def test_fail_scm():
    """ Test that both ADF and DFTB returns ``None`` if a computation fails"""
    # 5 membered ring from which ozone will dissociate
    mol = Molecule("test/test_files/ethylene.xyz")

    # Some dftb specific settings
    dftb_set = Settings()
    dftb_set.specific.dftb.dftb.scc

    # Calculate the DFTB hessian
    opt_dftb = dftb(templates.geometry.overlay(dftb_set), mol,
                    job_name="failed_DFTB")
    # fail_adf = adf(None, opt_dftb.molecule, job_name="fail_adf")
    # result = run(fail_adf.molecule)
    result = run(opt_dftb.molecule)
    print(result)
    assert isNone(result)


def test_fail_dirac():
    """ Dirac package should return ``None`` if it fails """
    mol = Molecule("test/test_files/h2.xyz")
    s = Settings()
    job = dirac(s, mol, job_name="fail_dirac")
    result = run(job.energy)

    assert isNone(result)


def test_fail_gamess():
    """
    Gamess should return ``None`` if a calculation fails.
    """
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

    result = run(methanol_geometry.molecule)

    assert isNone(result)


def test_fail_orca():
    """ Orca package should returns ``None`` if the computation fails"""
    methanol = Molecule('test/test_files/methanol.xyz')

    s = Settings()
    s.specific.orca.main = "RKS The_Cow_Functional SVP Opt TightSCF SmallPrint"

    opt = orca(s, methanol, job_name='fail_orca')
    result = run(opt.molecule)

    assert isNone(result)


if __name__ == "__main__":
    test_fail_gamess()
