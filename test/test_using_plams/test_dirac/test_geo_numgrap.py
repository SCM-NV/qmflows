from noodles import gather
from scm.plams import Molecule
from qmflows import (dirac, run, Settings)
import pytest


@pytest.mark.xfail
@pytest.mark.slow
def test_numgrad():
    """
    Test Dirac numerical opt.
    """
    # read geometry
    h2 = Molecule('test/test_files/h2.xyz')

    # create a dirac job
    job = create_job("h2_opt", h2)

    properties = [job.energy, job.molecule]
    rs = run(gather(*properties))

    print("Energy:", rs[0])
    print("Molecule: ", rs[1])


def create_job(name, mol):
    """
    Create a minimal optimization Job
    """
    s = Settings()

    s.specific.dirac.dirac['WAVE FUNCTION']
    s.specific.dirac.dirac.OPTIMIZE
    s.specific.dirac.dirac.OPTIMIZE["_en"] = True
    s.specific.dirac.dirac.OPTIMIZE.NUMGRA
    s.specific.dirac.HAMILTONIAN["LEVY-LEBLOND"]
    s.specific.dirac.HAMILTONIAN.DFT = "LDA"

    s.specific.dirac.GRID.RADINT = "1.0D-9"
    s.specific.dirac.GRID.ANGINT = 15

    s.specific.dirac.INTEGRALS.READIN["UNCONT"]

    s.specific.dirac['WAVE FUNCTION']["scf"]

    s.specific.dirac.molecule.basis.default = "cc-pVDZ"

    s.specific.dirac.molecule.coordinates.units = "AU"

    return dirac(s, mol, job_name=name)


if __name__ == "__main__":
    test_numgrad()
