"""Mock utilities."""

import pytest
import numpy as np
from assertionlib import assertion
from noodles import run_single
from pytest_mock import MockFixture
from typing import Any, List, Optional

from qmflows.test_utils import HAS_RDKIT

if HAS_RDKIT:
    from qmflows.components.operations import select_max, select_min


def generate_mocked_results(
        mocker: MockFixture, target: str, instances: int=10,
        expected: Optional[int]=None) -> List[Any]:
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


@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_select_max_list(mocker: MockFixture):
    """Test select_max using mocked results."""
    results = generate_mocked_results(
        mocker, 'qmflows.packages.SCM.ADF_Result', expected=1e3)
    wf = select_max(results, prop='prop')
    xs = run_single(wf)

    assertion.eq(xs.prop, 1e3)


@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_select_min(mocker: MockFixture):
    """Test select_min using mocked results."""
    results = generate_mocked_results(
        mocker, 'qmflows.packages.SCM.DFTB_Result', expected=-1e3)
    wf = select_min(results, prop='prop')

    xs = run_single(wf)

    assertion.eq(xs.prop, -1e3)
