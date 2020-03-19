"""Test CP2K parser functions."""
from assertionlib import assertion

from qmflows.parsers.cp2KParser import parse_cp2k_warnings, readCp2KBasis
from qmflows.test_utils import PATH
from qmflows.warnings_qmflows import QMFlows_Warning, cp2k_warnings


def test_parse_cp2k_warnings():
    """Parse CP2K warnings."""
    OUTPUT_FILE = PATH / "output_cp2k" / "cp2k_job" / "cp2k_job.out"
    map_warns = parse_cp2k_warnings(OUTPUT_FILE, cp2k_warnings)
    assertion.truth(all((val == QMFlows_Warning)
                        for val in map_warns.values()))


def test_read_basis():
    """Test that the basis are read correctly."""
    BASIS_FILE = PATH / "BASIS_MOLOPT"

    for key, data in zip(*readCp2KBasis(BASIS_FILE)):
        # The formats contains a list
        assertion.len(key.basisFormat)
        # Atoms are either 1 or two characters
        assertion.le(len(key.atom), 2)
        # All basis are MOLOPT
        assertion.contains(key.basis, "MOLOPT")
        # There is a list of exponents and coefficinets
        assertion.len(data.exponents)
        assertion.len(data.coefficients[0])
