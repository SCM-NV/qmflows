from nose.plugins.attrib import attr


@attr('slow')
def test_calc_freqs():
    """
    Test conditional workflow for freq calculation.
    """
    local_env = {}
    exec(open('examples/Conditional_workflows/calc_freqs.py').read(), {}, local_env)
    ref_table = "pbe         1533.267  3676.165  3817.097\n" +\
                "b3lyp       1515.799  3670.390  3825.813\n" +\
                "blyp        1529.691  3655.573  3794.110\n"
    assert str(local_env['table']) == ref_table
