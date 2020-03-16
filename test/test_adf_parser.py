"""Check methods to read adf output."""

from qmflows.test_utils import PATH
from qmflows.parsers.adf_parser import kfreader
from qmflows.parsers.generic_parsers import awk_file, extract_line_value
import numpy as np

path_t21 = PATH / "output_adf" / "ADFjob" / "ADFjob.t21"
path_output = PATH / "output_adf" / "ADFjob" / "ADFjob.out"


def test_adf_parser():
    """Test the adf output readers."""
    energy = kfreader(path_t21, section="Energy", prop="Bond Energy")
    charge = kfreader(path_t21, section="Properties",
                      prop="AtomCharge Mulliken")
    assert np.isclose(np.sum(charge), 0)
    assert np.isfinite(energy)


def test_adf_awk():
    """Test the awk parser for ADF."""
    script = r"/Total Used/ {print $9}"
    filename = path_output.as_posix()
    results = np.array(awk_file(filename, script=script))
    assert np.isfinite(np.sum(results))


def test_line_extraction():
    """Test the line extraction function for ADF."""
    homo = extract_line_value(path_output, pattern="HOMO", pos=4)
    assert np.isfinite(homo)
