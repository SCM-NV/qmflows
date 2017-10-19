
from scm.plams import (Atom, Molecule)
from qmflows import (dirac, Settings)


def main():

    pass
    # name = "BaF_req"
    # baf = create_molecule(name)

    # job = create_job(name, baf)

    # run(job)


def create_job(name, mol, basis='dyall.cv2z'):

    s = Settings()

    s.specific.dirac.dirac['WAVE FUNCTION']
    s.specific.dirac.dirac["4INDEX"]

    s.specific.dirac.HAMILTONIAN.X2Cmmf
    s.specific.dirac.HAMILTONIAN.GAUNT

    s.specific.dirac.INTEGRALS.READIN["UNCONT"]

    s.specific.dirac['WAVE FUNCTION']["relccsd"]
    s.specific.dirac['WAVE FUNCTION']["scf"]
    s.specific.dirac["WAVE FUNCTION"]["SCF"]["_en"] = True
    s.specific.dirac["WAVE FUNCTION"]["SCF"]["CLOSED SHELL"] = 64
    s.specific.dirac["WAVE FUNCTION"]["SCF"]["MAXITR"] = 50
    s.specific.dirac["WAVE FUNCTION"]["SCF"]["EVCCNV"] = "1.0D-9 5.0D-8"

    s.specific.dirac.moltra.active  = 'energy -3.0 30.0 0.01'

    s.specific.dirac.molecule.basis.default = basis
    s.specific.dirac.molecule.basis.special = 'F BASIS cc-pVDZ'

    s.specific.dirac.relccsd.FOCKSPACE
    s.specific.dirac.relccsd.CCSORT.NORECM
    s.specific.dirac.relccsd.CCFSPC.DOEA
    s.specific.dirac.relccsd.CCFSPC.NACTP = '6 6'

    return  dirac(s, mol, job_name=name)


def create_molecule(name, r=2.25):
    mol = Molecule()
    mol.add_atom(Atom(symbol='Ba', coords=(0, 0, 0)))
    mol.add_atom(Atom(symbol='F', coords=(0, 0, r)))

    return mol


if __name__ == "__main__":
    main()
