"""Test CP2K components."""
from qmflows.parsers import cp2KParser
from qmflows.test_utils import PATH


def test_read_cp2k_coefficients():
    """Test that the CP2K output coefficients are read properly."""
    path = PATH / "output_cp2k"
    mo_file = path / "mo_coeffs.log"

    data = cp2KParser.read_cp2k_coefficients(mo_file.as_posix(), path)

    assert data.eigenVals.ndim == 1 and data.coeffs.ndim == 2
