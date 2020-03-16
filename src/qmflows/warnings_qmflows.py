"""A module with warnings used throughout QMFlows."""

from .type_hints import WarnDict

__all__ = [
    'QMFlows_Warning', 'Generic_Key_Warning',
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


cp2k_warnings: WarnDict = {
    'SCF run NOT converged': SCF_Convergence_Warning
}
