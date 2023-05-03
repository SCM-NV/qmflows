"""A module with warnings used throughout QMFlows."""

from __future__ import annotations

from math import isclose
from collections.abc import Iterable

from pyparsing import ZeroOrMore, Suppress, SkipTo

from .type_hints import WarnMap
from .common import ParseWarning

__all__ = [
    'QMFlows_Warning', 'Key_Warning', 'Assertion_Warning',
    'SCF_Convergence_Warning', 'Geometry_Convergence_Warning',
    'Parameter_Warning', 'Charge_Warning', 'Orbital_Warning',
    'QMFlowsDeprecationWarning',
    'cp2k_warnings'
]


class QMFlows_Warning(Warning):
    """Generic :exc:`Warning` class for QMFlows."""


class Assertion_Warning(QMFlows_Warning):
    """Warning class for ``assertion``-related warnings."""


class Key_Warning(QMFlows_Warning):
    """Warning class for issues with :class:`~collections.abc.Mapping` keys."""


class SCF_Convergence_Warning(QMFlows_Warning):
    """Warning class for SCF-convergence issues."""


class Geometry_Convergence_Warning(QMFlows_Warning):
    """Warning class for geometry-convergence issues."""


class Parameter_Warning(QMFlows_Warning):
    """Warning class for classical forcefield parameters."""


class Charge_Warning(Parameter_Warning):
    """Warning class for charges in classical forcefields."""


class Orbital_Warning(QMFlows_Warning):
    """Warning class for orbital-related issues."""


class QMFlowsDeprecationWarning(DeprecationWarning, QMFlows_Warning):
    """Warning class for deprecations."""


def _eval_charge(msg: str, tolerance: float = 0.1) -> None | str:
    """Check of the total molecular charge is integer within a given *tolerance*."""
    charge = float(msg.rsplit(maxsplit=1)[1])
    charge_int = round(charge)

    condition = isclose(charge, charge_int, rel_tol=0, abs_tol=tolerance)
    return None if condition else msg.rstrip()


def _eval_param(
    msg: str,
    skip: Iterable[str] = ('Urey-Bradley', 'Out of plane bend'),
) -> None | str:
    """Return missing forcefield warnings in *msg* except for all terms in *skip*."""
    for i in skip:
        if i in msg:
            return None
    return msg.rstrip()


cp2k_warnings: WarnMap = {
    'SCF run NOT converged': ParseWarning(
        warn_type=SCF_Convergence_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("*** WARNING")) + SkipTo('\n\n'))
    ),

    'Missing': ParseWarning(
        warn_type=Parameter_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("FORCEFIELD| Missing")) + SkipTo('\n')),
        func=_eval_param
    ),

    'Total Charge of the Classical System:': ParseWarning(
        warn_type=Charge_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("CHARGE_INFO| Total Charge of "
                                          "the Classical System:")) + SkipTo('\n')),
        func=_eval_charge
    ),
}
