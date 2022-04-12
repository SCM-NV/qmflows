"""Deprecated main-namespace aliases for :mod:`qmflows.examples` functions."""

import functools
import warnings

from . import (
    example_freqs,
    example_H2O2_TS,
    example_generic_constraints,
    example_partial_geometry_opt,
)
from ..warnings_qmflows import QMFlowsDeprecationWarning

__all__ = [
    "_example_freqs",
    "_example_H2O2_TS",
    "_example_generic_constraints",
    "_example_partial_geometry_opt",
]


@functools.wraps(example_freqs)
def _example_freqs(*args, n_processes=1, **kwargs):
    """Deprecated alias for :func:`~qmflows.examples.example_freqs`."""
    warnings.warn(
        "`qmflows.example_freqs` is a deprecated alias for "
        "`qmflows.examples.example_freqs`",
        QMFlowsDeprecationWarning, stacklevel=2,
    )
    return example_freqs(*args, n_processes=n_processes, **kwargs)


@functools.wraps(example_H2O2_TS)
def _example_H2O2_TS(*args, **kwargs):
    """Deprecated alias for :func:`~qmflows.examples.example_H2O2_TS`."""
    warnings.warn(
        "`qmflows.example_H2O2_TS` is a deprecated alias for "
        "`qmflows.examples.example_H2O2_TS`",
        QMFlowsDeprecationWarning, stacklevel=2,
    )
    return example_H2O2_TS(*args, **kwargs)


@functools.wraps(example_generic_constraints)
def _example_generic_constraints(*args, n_processes=1, **kwargs):
    """Deprecated alias for :func:`~qmflows.examples.example_generic_constraints`."""
    warnings.warn(
        "`qmflows.example_generic_constraints` is a deprecated alias for "
        "`qmflows.examples.example_generic_constraints`",
        QMFlowsDeprecationWarning, stacklevel=2,
    )
    return example_generic_constraints(*args, n_processes=n_processes, **kwargs)


@functools.wraps(example_partial_geometry_opt)
def _example_partial_geometry_opt(*args, n_processes=1, **kwargs):
    """Deprecated alias for :func:`~qmflows.examples.example_partial_geometry_opt`."""
    warnings.warn(
        "`qmflows.example_partial_geometry_opt` is a deprecated alias for "
        "`qmflows.examples.example_partial_geometry_opt`",
        QMFlowsDeprecationWarning, stacklevel=2,
    )
    return example_partial_geometry_opt(*args, n_processes=n_processes, **kwargs)
