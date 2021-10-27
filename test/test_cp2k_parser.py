"""Test CP2K parser functions."""
from assertionlib import assertion
from distutils.spawn import find_executable

import pytest
from qmflows.parsers.cp2KParser import parse_cp2k_warnings, readCp2KBasis, get_cp2k_version_run
from qmflows.test_utils import PATH
from qmflows.warnings_qmflows import QMFlows_Warning, cp2k_warnings


def cp2k_available() -> bool:
    """Check if cp2k is installed."""
    path = find_executable("cp2k.popt")
    return path is not None


HAS_CP2K = cp2k_available()


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



@pytest.mark.skipif(not HAS_CP2K, reason="CP2K is not install or not loaded")
def test_get_cp2k_version() -> None:
    out = get_cp2k_version_run(PATH / "cp2k_freq.run")
    assertion.truth(out)
