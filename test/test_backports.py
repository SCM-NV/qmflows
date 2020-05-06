"""Tests for :mod:`qmflows.backports`."""

import sys
from typing import Type, Union

from assertionlib import assertion

from qmflows.backports import (Literal, Final, nullcontext,
                               _LiteralBackup, _FinalBackup, _NullContextBackup)


def test_literal():
    """Tests for :data:`Literal` and :data:`_LiteralBackup`."""
    assertion.eq(Literal.__module__, 'qmflows.backports')

    literal = _LiteralBackup()
    assertion.eq(literal[1], int)
    assertion.eq(literal[int], Type[int])
    assertion.eq(literal[1, 2], Union[int])
    assertion.eq(literal[1, 2.0, int, float], Union[int, float, Type[int], Type[float]])


def test_final():
    """Tests for :data:`Final` and :data:`_FinalBackup`."""
    assertion.eq(Final.__module__, 'qmflows.backports')

    final = _FinalBackup()
    assertion.eq(final[int], int)

    assertion.assert_(final.__getitem__, 1, exception=TypeError)
    assertion.assert_(final.__getitem__, (int, float), exception=TypeError)


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
