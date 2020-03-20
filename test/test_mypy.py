"""Tests using mypy."""

import warnings
import subprocess
from os.path import dirname, join
from typing import Tuple
from pathlib import Path

from qmflows.test_utils import PATH, Assertion_Warning
from qmflows.backports import Literal

_ROOT = Path(__file__).parts[:-2]
PACKAGE = join(*_ROOT, 'src', 'qmflows')
INI = str(PATH / 'mypy.ini')

#: Type annotation for the 'action' keyword.
Action = Literal['raise', 'warn', 'ignore']


def test_mypy(action: Action = 'warn') -> None:
    """Test using mypy."""
    command = f"mypy {PACKAGE!r} --config-file {INI!r}"
    out = subprocess.run(command, capture_output=True, shell=True)

    stdout = out.stdout.decode()

    if action == 'warn' and out.returncode != 0:
        warnings.warn(stdout, category=Assertion_Warning)

    elif action == 'raise':
        try:
            assert out.returncode == 0, stdout
        except AssertionError as ex:
            msg = stdout.rsplit('\n', maxsplit=2)[1]
            raise AssertionError(msg) from ex

    else:
        pass

