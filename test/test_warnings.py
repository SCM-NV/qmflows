"""Tests for :mod:`qmflows.warnings_qmflows`."""

from assertionlib import assertion

from qmflows.warnings_qmflows import (
    QMFlows_Warning, Key_Warning, SCF_Convergence_Warning, Assertion_Warning,
    Geometry_Convergence_Warning, Parameter_Warning, Charge_Warning,
    _eval_charge, _eval_param
)


def test_warnings() -> None:
    """Tests for all QMFlows :exc:`Warning` types."""
    assertion.issubclass(QMFlows_Warning, Warning)
    assertion.issubclass(Assertion_Warning, QMFlows_Warning)
    assertion.issubclass(Key_Warning, QMFlows_Warning)
    assertion.issubclass(SCF_Convergence_Warning, QMFlows_Warning)
    assertion.issubclass(Geometry_Convergence_Warning, QMFlows_Warning)
    assertion.issubclass(Parameter_Warning, QMFlows_Warning)
    assertion.issubclass(Charge_Warning, Parameter_Warning)


def test_eval_charge() -> None:
    """Tests for :func:`_eval_charge`."""
    msg1 = 'bob  1.2    '
    msg2 = 'bob  1.02   '

    out1 = _eval_charge(msg1)
    out2 = _eval_charge(msg2)

    assertion.eq(out1, msg1.rstrip())
    assertion.is_(out2, None)


def test_eval_param() -> None:
    """Tests for :func:`_eval_param`."""
    msg1 = 'Angles          '
    msg2 = 'Urey-Bradley    '

    out1 = _eval_param(msg1)
    out2 = _eval_param(msg2)

    assertion.eq(out1, msg1.rstrip())
    assertion.is_(out2, None)
