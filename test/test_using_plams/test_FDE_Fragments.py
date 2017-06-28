from nose.plugins.attrib import attr
from qmworks.examples import (
    example_ADF3FDE_Cystine, example_ADF3FDE_Dialanine, example_FDE_fragments)
import numpy as np


@attr('slow')
def test_ADF3FDE_Dialanine():
    """
    Test MFCC partitioning of dialanine
    """
    expected_dipoles = np.array(([2.2836628592, 1.6858674039, 0.0719124986],
                                 [2.28429496, 1.68676037, 0.07201473]))
    assert_dipoles(expected_dipoles, example_ADF3FDE_Dialanine())


@attr('slow')
def test_ADF3FDE_Cystine():
    """
    Test MFCC partitioning of cystine
    """
    expected_dipoles = np.array(
        ([-0.3730226317861707, -0.43888610949488793, 0.14005294054732076],
         [-0.35257865, -0.46682244, 0.17783698],
         [-0.36767441, -0.46617708, 0.16977581]))

    assert_dipoles(expected_dipoles, example_ADF3FDE_Cystine())


@attr('slow')
def test_FDE_Fragments():
    """
    Test FDE Fragments
    """
    expected_dipoles = np.array([7.12156395e-01, 1.22645358e-08, 2.50112616e-17])
    test_dipoles = example_FDE_fragments()
    print(test_dipoles)
    assert_dipoles(expected_dipoles, test_dipoles)


def assert_dipoles(expected, test, rtol=1e-3):
    """
    Check that the expected dipoles is close to the example.
    """
    assert np.allclose(expected, test, rtol=rtol)
