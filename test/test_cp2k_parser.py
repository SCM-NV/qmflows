"""Test CP2K parser functions."""
from assertionlib import assertion

from qmflows.parsers.cp2KParser import parse_cp2k_warnings
from qmflows.test_utils import PATH
from qmflows.warnings_qmflows import QMFlows_Warning, cp2k_warnings

OUTPUT_FILE = PATH / "output_cp2k" / "cp2k_job" / "cp2k_job.out"


def test_parse_cp2k_warnings():
    """Parse CP2K warnings."""
    map_warns = parse_cp2k_warnings(OUTPUT_FILE, cp2k_warnings)
    assertion.truth(all((val == QMFlows_Warning) for val in map_warns.values()))
