from nose.plugins.attrib import attr


@attr('slow')
def test_calc_freqs():
    """
    Test conditional workflow for freq calculation.
    """
    local_env = {}
    exec(open('examples/Conditional_workflows/calc_freqs.py').read(), {}, local_env)
    ref_table = "pbe         1532.932  3677.915  3818.842\n" +\
                "b3lyp       1515.407  3672.076  3827.516\n" +\
                "blyp        1529.359  3657.323  3795.853\n"
    assert str(local_env['table']) == ref_table
