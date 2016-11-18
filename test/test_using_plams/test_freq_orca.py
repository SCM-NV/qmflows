
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks import (orca, run, Settings)


@attr('slow')
def test_freq_ethylene():
    """
    Run a methanol optimization and retrieve the optimized geom.
    """
    ethylene = Molecule('test/test_files/ethylene.xyz')

    s = Settings()
    s.specific.orca.main = "freq"
    s.specific.orca.basis.basis = 'sto_sz'
    s.specific.orca.method.functional = 'lda'
    s.specific.orca.method.method = 'dft'

    freq = orca(s, ethylene)

    rs = run(freq.frequencies)
    print(rs)
