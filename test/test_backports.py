"""Tests for :mod:`qmflows.backports`."""

import sys
from typing import Type, Any

from assertionlib import assertion

from qmflows.backports import (
    Literal, Final, nullcontext, _NullContextBackup
)


def test_literal():
    """Tests for :data:`Literal` and :data:`_LiteralBackup`."""
    assertion.is_(Literal[1], Any)


def test_final():
    """Tests for :data:`Final` and :data:`_FinalBackup`."""
    assertion.is_(Final[int], Any)


def test_nullcontext():
    """Tests for :class:`nullcontext` and :data:`_NullContextBackup`."""
    if sys.version_info.minor >= 7:
        assertion.eq(nullcontext.__module__, 'contextlib')
    else:
        assertion.eq(nullcontext.__module__, 'qmflows.backports')

    cm = _NullContextBackup(True)
    with cm as bool_:
        assertion.is_(bool_, True)
    with cm as bool_:
        assertion.is_(bool_, True)
