"""Test CP2K parser functions."""

import os
import shutil
from pathlib import Path

import pytest
from assertionlib import assertion

from qmflows.parsers.cp2KParser import parse_cp2k_warnings, readCp2KBasis, get_cp2k_version_run
from qmflows.test_utils import PATH, requires_cp2k
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


@requires_cp2k
class TestGetCP2KVersion:
    @classmethod
    def setup_class(cls) -> None:
        """Add ``cp2k.popt`` to a space-containing directory and add it to ``$PATH``."""
        cp2k_popt = shutil.which("cp2k.popt")
        assert cp2k_popt is not None

        os.mkdir(PATH / "test dir")
        shutil.copy2(cp2k_popt, PATH / "test dir" / "cp2k.popt")

    @classmethod
    def teardown_class(cls) -> None:
        """Restore the old ``$PATH`` environment variable."""
        shutil.rmtree(PATH / "test dir")

    PASS_DICT = {
        "mpirun": "bin/mpirun cp2k.popt -i file.in -o file.out",
        "srun": "srun cp2k.popt -i file.in -o file.out",
        "none": "cp2k.popt -i file.in -o file.out",
        "space": " cp2k.popt -i file.in -o file.out",
        "single_quote": f"'{PATH}/test dir/cp2k.popt' -i file.in -o file.out",
        "double_quote": f'"{PATH}/test dir/cp2k.popt" -i file.in -o file.out',
    }

    @pytest.mark.parametrize("name,txt", PASS_DICT.items(), ids=PASS_DICT)
    def test_pass(self, name: str, txt: str, tmp_path: Path) -> None:
        with open(tmp_path / f"{name}.run", "w") as f:
            f.write(txt)
        out = get_cp2k_version_run(tmp_path / f"{name}.run")
        assertion.truth(out)

    RAISE_DICT = {
        "invalid_run": "cp2k.popt -j file.in -o file.out",
        "invalid_exec": "cp2k.bob -i file.in -o file.out",
        "invalid_exec_version": "python -i file.in -o file.out",
    }

    @pytest.mark.parametrize("name,txt", RAISE_DICT.items(), ids=RAISE_DICT)
    def test_raise(self, name: str, txt: str, tmp_path: Path) -> None:
        with open(tmp_path / f"{name}.run", "w") as f:
            f.write(txt)
        with pytest.raises(ValueError):
            get_cp2k_version_run(tmp_path / f"{name}.run")
