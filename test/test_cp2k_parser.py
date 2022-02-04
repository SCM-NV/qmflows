"""Test CP2K parser functions."""

import os
import shutil
from pathlib import Path

import numpy as np
import h5py
import pytest
from assertionlib import assertion

from qmflows.parsers.cp2KParser import parse_cp2k_warnings, readCp2KBasis, get_cp2k_version_run
from qmflows.test_utils import PATH, requires_cp2k
from qmflows.warnings_qmflows import QMFlows_Warning, cp2k_warnings
from qmflows.common import AtomBasisKey


def test_parse_cp2k_warnings():
    """Parse CP2K warnings."""
    OUTPUT_FILE = PATH / "output_cp2k" / "cp2k_job" / "cp2k_job.out"
    map_warns = parse_cp2k_warnings(OUTPUT_FILE, cp2k_warnings)
    assertion.truth(all((val == QMFlows_Warning)
                        for val in map_warns.values()))


class TestReadBasis:
    """Test that the basis are read correctly."""

    @staticmethod
    def get_key(key_tup: AtomBasisKey) -> str:
        return os.path.join(
            key_tup.atom,
            key_tup.basis,
            "-".join(str(i) for i in key_tup.basisFormat),
        )

    def test_pass(self) -> None:
        basis_file = PATH / "BASIS_MOLOPT"
        keys, values = readCp2KBasis(basis_file)

        with h5py.File(PATH / "basis.hdf5", "r") as f:
            for key_tup, value_tup in zip(keys, values):
                key = self.get_key(key_tup)
                exponents = f[key]["exponents"][...].T
                coefficients = f[key]["coefficients"][...].T
                np.testing.assert_allclose(value_tup.exponents, exponents, err_msg=key)
                np.testing.assert_allclose(value_tup.coefficients, coefficients, err_msg=key)

    PARAMS = {
        "BASIS_INVALID_N_EXP": 11,
        "BASIS_INVALID_FMT": 3,
        "BASIS_NOTIMPLEMENTED": 2,
    }

    @pytest.mark.parametrize("filename,lineno", PARAMS.items(), ids=PARAMS)
    def test_raise(self, filename: str, lineno: int) -> None:
        pattern = r'Failed to parse the basis set ".+" on line {}'.format(lineno)
        with pytest.raises(ValueError, match=pattern):
            readCp2KBasis(PATH / filename)


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
