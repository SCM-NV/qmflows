"""Tests for :mod:`qmflows.utils`."""

import textwrap
from io import TextIOBase, StringIO
from pathlib import Path
from contextlib import AbstractContextManager

import pytest
from assertionlib import assertion
from scm.plams import init, finish

from qmflows import Settings
from qmflows.packages.cp2k_package import CP2K_Result
from qmflows.utils import to_runtime_error, file_to_context, init_restart, InitRestart
from qmflows.test_utils import PATH_MOLECULES, PATH, validate_status

PSF_STR: str = """
PSF EXT

    10 !NATOM
     1 MOL1     1        LIG      C        C331   -0.272182       12.010600        0
     2 MOL1     1        LIG      C        C321   -0.282182       12.010600        0
     3 MOL1     1        LIG      C        C2O3    0.134065       12.010600        0
     4 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
     5 MOL1     1        LIG      O        O2D2   -0.210848       15.999400        0
     6 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
     7 MOL1     1        LIG      H        HGA2    0.087818        1.007980        0
     8 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
     9 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
    10 MOL1     1        LIG      H        HGA3    0.087818        1.007980        0
"""


def _test_function(settings, key, value, mol):
    raise Exception('test')


def test_to_runtime_error() -> None:
    """Tests for :func:`to_runtime_error`."""
    args = (None, '_test_function', None, None)
    func = to_runtime_error(_test_function)
    assertion.assert_(func, *args, exception=RuntimeError)


def test_file_to_context() -> None:
    """Tests for :func:`file_to_context`."""
    path_like = PATH_MOLECULES / 'mol.psf'
    file_like = StringIO(PSF_STR)

    cm1 = file_to_context(path_like)
    cm2 = file_to_context(file_like)
    assertion.isinstance(cm1, AbstractContextManager)
    assertion.isinstance(cm2, AbstractContextManager)
    assertion.isinstance(cm1.__enter__(), TextIOBase)
    assertion.isinstance(cm2.__enter__(), TextIOBase)

    sequence = range(10)
    assertion.assert_(file_to_context, sequence, require_iterator=False)
    assertion.assert_(file_to_context, sequence, require_iterator=True, exception=TypeError)

    assertion.assert_(file_to_context, None, exception=TypeError)
    assertion.assert_(file_to_context, [1, 2], exception=TypeError)
    assertion.assert_(file_to_context, 5.0, exception=TypeError)


def test_restart_init(tmp_path: Path) -> None:
    """Tests for :func:`restart_init` and :class:`RestartInit`."""
    workdir = tmp_path / 'test_restart_init'

    init(path=tmp_path, folder="test_restart_init")
    finish()
    assertion.isdir(workdir)

    init_restart(path=tmp_path, folder="test_restart_init")
    assertion.isdir(workdir)
    assertion.isdir(f'{workdir}.002', invert=True)
    finish()

    with InitRestart(path=tmp_path, folder="test_restart_init"):
        assertion.isdir(workdir)
        assertion.isdir(f'{workdir}.002', invert=True)


def test_validate_status() -> None:
    root = PATH / "output_cp2k" / "cp2k_opt"
    result = CP2K_Result(
        settings=Settings(),
        molecule=None,
        job_name="cp2k_opt",
        dill_path=root / "cp2k_opt.dill",
        plams_dir=root,
        work_dir=root,
        status="failed",
    )

    with pytest.raises(AssertionError) as rec:
        validate_status(result)
    msg = str(rec.value)

    with open(root / "cp2k_opt.out", "r") as f:
        out_ref = textwrap.indent("".join(f.readlines()[-100:]), 4 * " ")
    with open(root / "cp2k_opt.err", "r") as f:
        err_ref = textwrap.indent("".join(f.readlines()[-100:]), 4 * " ")

    assertion.contains(msg, out_ref)
    assertion.contains(msg, err_ref)
    assertion.contains(msg, f"Unexpected {result.job_name} status: {result.status!r}")
