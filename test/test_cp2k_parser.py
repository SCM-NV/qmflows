"""Test CP2K components."""
from pathlib import Path
from qmflows.parsers import cp2KParser


def test_read_cp2k_coefficients():
    """Test that the CP2K output coefficients are read properly."""
    root = Path("test")
    path = root / "test_files" / "output_cp2k"
    mo_file = path / "mo_coeffs.log"

    data = cp2KParser.read_cp2k_coefficients(mo_file.as_posix(), path)

    assert data.eigenVals.ndim == 1 and data.coeffs.ndim == 2
