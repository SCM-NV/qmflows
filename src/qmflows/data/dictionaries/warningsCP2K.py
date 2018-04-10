__all__ = ['cp2k_warnings']

from qmflows.warnings_qmflows import SCF_Convergence_Warning

cp2k_warnings = {
    'SCF run NOT converged': SCF_Convergence_Warning
}
