
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks.packages.orca import orca
from qmworks.packages import (run, Settings)


@attr('slow')
def test_freq_ethylene():
    """
    Run a methanol optimization and retrieve the optimized geom.
    """
    ethylene = Molecule('test/test_files/ethylene.xyz')

    s = Settings()
    s.specific.orca.main = " B3LYP SVP NumFreq SmallPrint"

    freq = orca(s, ethylene)

    mol = run(freq.molecule)
    print(mol)
    
    assert False
