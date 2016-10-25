from nose.plugins.attrib import attr
from test import exec_example

@attr('slow')
def test_ADF3FDE_Dialanine():
    """
    Test MFCC partitioning of dialanine
    """
    local_env = exec_example('examples/FDE_Fragments', 'ADF3FDE_Dialanine.py')
    for x,y in zip(local_env['supermol_dipole'], local_env['mfcc_dipole']):
        assert abs(x-y) < 0.01

@attr('slow')
def test_ADF3FDE_Cystine():
    """
    Test MFCC partitioning of cystine
    """
    local_env = exec_example('examples/FDE_Fragments', 'ADF3FDE_Cystine.py')
    for x,y in zip(local_env['supermol_dipole'], local_env['mfcc_dipole']):
        assert abs(x-y) < 0.01

@attr('slow')
def test_FDE_Fragments():
    """
    Test FDE Fragments
    """
    local_env = exec_example('examples/FDE_Fragments', 'FDE_Fragments.py')
    dipole_vect = local_env['dipole_fde']
    assert round(dipole_vect[0] - 0.71, 2) == 0
    assert round(dipole_vect[1], 2) == 0
    assert round(dipole_vect[2], 2) == 0
