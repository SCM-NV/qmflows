"""Check conditional workflows funcionality."""
import numpy as np
import pytest

from qmflows.test_utils import HAS_RDKIT, requires_orca, requires_adf

if HAS_RDKIT:
    from qmflows.examples import example_freqs


@pytest.mark.slow
@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
@requires_orca
@requires_adf
def test_calc_freqs():
    """Test conditional workflow for freq calculation."""
    test = example_freqs()
    expected = np.array(
        [[1533.26703326, 3676.16470838, 3817.0971787],
         [1515.798647, 3670.390391, 3825.813363],
         [1529.69058854, 3655.57330642, 3794.10972377]])

    assert np.allclose(test, expected)
