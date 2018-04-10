from nose.plugins.attrib import attr
from qmflows.parsers.gamess_parser import (parse_frequencies, parse_hessian)
import numpy as np

file_name = 'test/test_files/CNH_hess.dat'


@attr('fast')
def test_hess_gamess():
    """
    Test if the hessian is read properly.
    """
    hess = parse_hessian(file_name)
    hess = hess.reshape(9, 9)

    expected = [1.10665699e+00, -2.42593102e-02, -2.03554235e-08,
                1.02817225e+00, 1.52947465e-01, -1.60761763e-07,
                -4.86900429e-02, 3.41009491e-01, 1.51937485e-09]

    assert abs(np.sum(hess.diagonal() - expected)) < 1.0e-7


@attr('fast')
def test_frequencies_gamess():
    """
    Read if the frequencies are read from the *.dat file
    """
    freqs_modes = parse_frequencies(file_name)
    expected_freqs = [1.24863666e+03, 2.31217000, 1.46703000, 6.99020000e-01,
                      1.00000000e-03, 2.51780000e-01, 2.34140000e+00, 2.10562417e+03,
                      3.07137245e+03]
    assert abs(np.sum(freqs_modes[0] - expected_freqs)) < 1e-7
