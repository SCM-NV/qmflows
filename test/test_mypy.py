"""Tests using mypy."""

import warnings
import subprocess
from os.path import join, abspath
from pathlib import Path

import pytest
from typing_extensions import Literal
from qmflows.test_utils import PATH, Assertion_Warning

try:
    import mypy
except ImportError:
    HAS_MYPY = False
else:
    HAS_MYPY = True

_ROOT = Path(__file__).parts[:-2]
PACKAGE = join(*_ROOT, 'src', 'qmflows')
INI = abspath(PATH / 'mypy.ini')

Action = Literal['raise', 'warn', 'ignore']  #: Type annotation for the 'action' keyword.
ACTION = frozenset(['raise', 'warn', 'ignore'])


@pytest.mark.skipif(not HAS_MYPY, reason="Requires MyPy")
def test_mypy(action: Action = 'warn') -> None:
    """Test using mypy."""
    if action not in ACTION:
        raise ValueError(f"Invalid value for the 'action' parameter: {action!r}")

    command = f"mypy {PACKAGE!r} --config-file {INI!r}"
    out = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout = out.stdout.decode()
    stderr = out.stderr.decode()
    assert not stderr, stderr

    if action == 'warn' and out.returncode != 0:
        warnings.warn(stdout, category=Assertion_Warning)

    elif action == 'raise':
        try:
            assert out.returncode == 0, stdout
        except AssertionError as ex:
            msg = stdout.rsplit('\n', maxsplit=2)[1]
            raise AssertionError(msg) from ex
