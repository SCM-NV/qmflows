from nose.plugins.attrib import attr
from test import exec_example


@attr('slow')
def test_partial_geometry_opt():
    """
    Test partial geometry optimization.
    """
    local_env = exec_example('examples/Constrained_and_TS_optimizations', 'partial_geometry_opt.py')
    assert str(local_env['geom1']) == str(local_env['geom2'])
