"""Tests using mypy."""

import subprocess
from os.path import dirname, join
from typing import Tuple
from pathlib import Path

from qmflows.test_utils import PATH

_ROOT = Path(__file__).parts[:-2]
PACKAGE = join(*_ROOT, 'src', 'qmflows')
INI = str(PATH / 'mypy.ini')


def test_mypy() -> None:
    """Test using mypy."""
    command = f"mypy {PACKAGE!r} --config-file {INI!r}"
    out = subprocess.run(command, capture_output=True, shell=True)

    stdout = out.stdout.decode()
    try:
        assert out.returncode == 0, stdout
    except AssertionError as ex:
        msg = stdout.rsplit('\n', maxsplit=2)[1]
        raise AssertionError(msg) from ex
