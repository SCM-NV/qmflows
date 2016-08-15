from noodles import gather
from nose.plugins.attrib import attr
from qmworks import (Settings, adf, run)
from plams import Molecule
import plams


@attr('slow')
def test_freq():
    """
    Do some constraint optimizations then launch a freq calc.
    """
    plams.init()

    job = geo_opt()
    rs = freq_calc(job.molecule)
    results = run(rs)
    plams.finish()

    assert False


def geo_opt():
    """
    create some optimizations jobs
    """
    mol = Molecule("test/test_files/ethene.xyz", "xyz")

    s = Settings()
    s.specific.adf.xc.GGA = "BP86"
    s.specific.adf.basis.type = "DZP"
    s.specific.adf.beckegrid.quality = "good"
    s.specific.adf.geometry.converge = "e=0.0001"

    s.specific.adf.geometry.optim = "delocalized"
    s.specific.adf.scf.iterations = "99"
    s.specific.adf.scf.converge = "0.0000001"
    s.specific.adf.constraints.dist = "1 2 1.4"

    return adf(s, mol)


def freq_calc(mol):
    """
    Read geometries from an optimization job and do freq calculations.
    """
    s = Settings()
    s.specific.adf.xc.GGA = "BP86"
    s.specific.adf.basis.type = "DZP"
    s.specific.adf.beckegrid.quality = "good"
    s.specific.adf.scf.iterations = "99"
    s.specific.adf.scf.converge = "0.0000001"
    s.specific.adf.analyticalfreq

    return adf(s, mol)
