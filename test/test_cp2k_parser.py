from qmflows.parsers import cp2KParser
import os


def test_read_cp2k_coefficients():
    """Test that the CP2K output coefficients are read properly."""
    path = "test/test_files/output_cp2k"
    mo_file = os.path.join(path, "mo_coeffs.log")

    data = cp2KParser.read_cp2k_coefficients(mo_file, path)

    assert data.eigenVals.ndim == 1 and data.coeffs.ndim == 2
