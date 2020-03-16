"""Check methods to read adf output."""

from qmflows.test_utils import PATH
from qmflows.parsers.adf_parser import kfreader
from qmflows.parsers.generic_parsers import awk_file
import numpy as np

path_t21 = PATH / "output_adf" / "ADFjob" / "ADFjob.t21"


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
    path_output = PATH / "output_adf" / "ADFjob"
    filename = (path_output / "ADFjob.out").as_posix()
    results = np.array(awk_file(filename, script=script))
    assert np.isfinite(np.sum(results))
