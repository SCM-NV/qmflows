

from qmworks.parsers.gamess_parser import parse_hess
import numpy as np

file_name = 'test/test_files/CNH_hess.dat'


def test_hess_gamess():
    """
    Test if the hessian is read properly.
    """
    hess = parse_hess(file_name)
    hess = hess.reshape(9, 9)

    expected = [1.10665699e+00, -2.42593102e-02, -2.03554235e-08,
                1.02817225e+00, 1.52947465e-01, -1.60761763e-07,
                -4.86900429e-02, 3.41009491e-01, 1.51937485e-09]

    assert np.sum(hess.diagonal() - expected) < 1.0e-7


def test_frequencies_gamess():
    """
    Read if the frequencies are read from the *.dat file
    """
