"""Mock utilities."""
import numpy as np
from assertionlib import assertion
from noodles import run_single
from pytest_mock import mocker

from qmflows.components.operations import select_max, select_min


def generate_mocked_results(mocker, target, instances=10, expected=None):
    """Generate a list of mocked results with property `prop`.

    One of the results is set to `expected`  and the rest are random values.
    """

    results = [mocker.patch(target) for _ in range(10)]

    # Set the last Results with the expected value
    results[9].prop = expected

    # Initialize the other results with random numbers
    rs = np.random.uniform(high=5.0, size=9)

    for (i,), x in np.ndenumerate(rs):
        setattr(results[i], 'prop', x)

    return results


def test_select_max_list(mocker):
    """Test select_max using mocked results."""
    results = generate_mocked_results(
        mocker, 'qmflows.packages.SCM.ADF_Result', expected=1e3)
    wf = select_max(results, prop='prop')
    xs = run_single(wf)

    assertion.eq(xs.prop, 1e3)


def test_select_min(mocker):
    """Test select_min using mocked results."""
    results = generate_mocked_results(
        mocker, 'qmflows.packages.SCM.DFTB_Result', expected=-1e3)
    wf = select_min(results, prop='prop')

    xs = run_single(wf)

    assertion.eq(xs.prop, -1e3)
