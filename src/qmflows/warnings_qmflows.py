"""A module with warnings used throughout QMFlows."""

from math import isclose, inf
from typing import Optional, Mapping, Collection

from pyparsing import ZeroOrMore, Suppress, SkipTo

from .type_hints import ParseWarning

__all__ = [
    'QMFlows_Warning', 'Key_Warning',
    'SCF_Convergence_Warning', 'Geometry_Convergence_Warning',
    'Parameter_Warning', 'Charge_Warning',
    'cp2k_warnings'
]


class QMFlows_Warning(Warning):
    """Generic :exc:`Warning` class for QMFlows."""


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


def _eval_charge(msg: str, tolerance: float = 0.1) -> Optional[str]:
    """Check of the total moleculair charge is integer within a given *tolerance*."""
    charge = float(msg.rsplit(' ', maxsplit=1)[1])
    charge_int = int(charge)

    condition = isclose(charge, charge_int, rel_tol=inf, abs_tol=tolerance)
    if condition:
        return None
    else:
        return msg.strip().rstrip()


def _eval_param(msg: str, skip: Collection[str] = ('Urey-Bradley', 'Out of plane bend')
                ) -> Optional[str]:
    """Return missing forcefield warnings in *msg* except for all terms in *skip*."""
    msg_gen = (i.strip().rstrip() for i in msg.splitlines() if 'FORCEFIELD| Missing' in i)
    ret = '\n'.join(i for i in msg_gen if all(j not in i for j in skip))
    return ret or None


def _return_msg(msg: str) -> str:
    """Return the passed *msg*."""
    return msg


cp2k_warnings: Mapping[str, ParseWarning] = {
    'SCF run NOT converged': ParseWarning(
        warn_type=SCF_Convergence_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("*** WARNING")) + SkipTo('\n\n')),
        func=_return_msg
    ),

    'Missing': ParseWarning(
        warn_type=Parameter_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("FORCEFIELD| WARNING:")) + SkipTo('Charge Checking')),
        func=_eval_param
    ),

    'Total Charge of the Classical System:': ParseWarning(
        warn_type=Charge_Warning,
        parser=ZeroOrMore(Suppress(SkipTo("CHARGE_INFO| Total Charge of the Classical System:")) + SkipTo('\n')),
        func=_eval_charge
    ),
}
