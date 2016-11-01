from nose.plugins.attrib import attr


@attr('slow')
def test_partial_geometry_opt():
    """
    Test partial geometry optimization.
    """
    local_env = {}
    exec(open('examples/Constrained_and_TS_optimizations/partial_geometry_opt.py').read(), {}, local_env)
    assert str(local_env['geom1']) == str(local_env['geom2'])

def test_13dipolar_cyclo_ts():
    """
    Test TS optimization of a 1,3-dipolar-cycloaddition in ORCA using init hessian from DFTB
    """
    local_env = {}
    exec(open('examples/Constrained_and_TS_optimizations/13DipolarCyclo_TS.py').read(), {}, local_env)
    assert local_env['optcycles'] < 5
    assert abs(local_env['d1'] - 2.294) < 0.01
    assert abs(local_env['d2'] - 2.294) < 0.01

def test_generic_constraints():
    """
    Test generic distance constraints on all packages
    """
    local_env = {}
    exec(open('examples/Constrained_and_TS_optimizations/generic_constraints.py').read(), {}, local_env)
    assert local_env['table'] == {'adf': {'1.1': -0.276479, '1.0': -0.284003, '1.2': -0.258743},
                                  'dftb': {'1.1': -4.747407, '1.0': -4.760192, '1.2': -4.732274},
                                  'orca': {'1.1': -99.430973, '1.0': -99.438656, '1.2': -99.415615}}
