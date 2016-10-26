from nose.plugins.attrib import attr
from test import exec_example


@attr('slow')
def test_partial_geometry_opt():
    """
    Test partial geometry optimization.
    """
    local_env = exec_example('examples/Constrained_and_TS_optimizations', 'partial_geometry_opt.py')
    assert str(local_env['geom1']) == str(local_env['geom2'])

def test_13dipolar_cyclo_ts():
    """
    Test TS optimization of a 1,3-dipolar-cycloaddition in ORCA using init hessian from DFTB
    """
    local_env = exec_example('examples/Constrained_and_TS_optimizations', '13DipolarCyclo_TS.py')
    assert local_env['optcycles'] < 5
    assert abs(local_env['d1'] - 2.294) < 0.01
    assert abs(local_env['d2'] - 2.294) < 0.01
