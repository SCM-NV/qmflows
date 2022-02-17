"""Test CP2K parser functions."""

import os
import shutil
from pathlib import Path

import numpy as np
import h5py
import pytest
from assertionlib import assertion

from qmflows.parsers.cp2KParser import parse_cp2k_warnings, readCp2KBasis
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
        pattern = r"Failed to parse the '.+' basis set on line {}".format(lineno)
        with pytest.raises(ValueError, match=pattern):
            readCp2KBasis(PATH / filename)
