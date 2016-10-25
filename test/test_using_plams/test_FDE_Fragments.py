from nose.plugins.attrib import attr
import os
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
    local_env = {}
    global_env = {}
    os.chdir('examples/FDE_Fragments')
    try:
        exec(open('ADF3FDE_Cystine.py').read(), global_env, local_env)
        os.chdir('../..')
    except:
        os.chdir('../..')
        raise RuntimeError()
    for x,y in zip(local_env['supermol_dipole'], local_env['mfcc_dipole']):
        assert abs(x-y) < 0.01

@attr('slow')
def test_FDE_Fragments():
    """
    Test FDE Fragments
    """
    local_env = {}
    global_env = {}
    os.chdir('examples/FDE_Fragments')
    try:
        exec(open('FDE_Fragments.py').read(), global_env, local_env)
        os.chdir('../..')
    except:
        os.chdir('../..')
        raise RuntimeError()
    dipole_vect = local_env['dipole_fde']
    assert round(dipole_vect[0] - 0.71, 2) == 0
    assert round(dipole_vect[1], 2) == 0
    assert round(dipole_vect[2], 2) == 0
