__all__ = ['SCF_Convergence_Warning', 'Geometry_Convergence_Warning', 'cp2k_warnings']


class SCF_Convergence_Warning(RuntimeWarning):
    pass


class Geometry_Convergence_Warning(RuntimeWarning):
    pass


cp2k_warnings = {
    'SCF run NOT converged': SCF_Convergence_Warning
}
