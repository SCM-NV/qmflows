"""Tests for :mod:`qmflows.warnings_qmflows`."""

import pytest
from assertionlib import assertion

from qmflows.warnings_qmflows import (
    QMFlows_Warning, Key_Warning, SCF_Convergence_Warning, Assertion_Warning,
    Geometry_Convergence_Warning, Parameter_Warning, Charge_Warning,
    Orbital_Warning, QMFlowsDeprecationWarning, _eval_charge, _eval_param,
)


@pytest.mark.parametrize("cls", [
    QMFlows_Warning,
    Key_Warning,
    Assertion_Warning,
    SCF_Convergence_Warning,
    Geometry_Convergence_Warning,
    Parameter_Warning,
    Charge_Warning,
    Orbital_Warning,
    QMFlowsDeprecationWarning,
])
def test_warnings(cls: "type[QMFlows_Warning]") -> None:
    """Tests for all QMFlows :exc:`Warning` types."""
    warning = cls()
    assertion.isinstance(warning, QMFlows_Warning)


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
