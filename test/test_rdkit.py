"""Tests related to the optional rdkit dependency."""

import sys
import importlib

import pytest
import qmflows
from qmflows.test_utils import HAS_RDKIT
from assertionlib import assertion

NAMES = frozenset({
    'components',
    'Angle',
    'Dihedral',
    'Distance',
    'find_first_job',
    'select_max',
    'select_min',
    'examples',
    'example_H2O2_TS',
    'example_freqs',
    'example_generic_constraints',
    'example_partial_geometry_opt',
})


@pytest.mark.parametrize("name", sorted(NAMES))
def test_sub_module(name: str) -> None:
    """Test :func:`getattr` operations on rdkit-requiring objects."""
    if HAS_RDKIT:
        assert getattr(qmflows, name)
    else:
        with pytest.raises(ImportError):
            getattr(qmflows, name)


@pytest.mark.skipif(sys.version_info < (3, 7), reason="requires python 3.7")
def test_requires_rdkit() -> None:
    """Test that ``NAMES`` and ``qmflows._REQUIRES_RDKIT`` are synced."""
    assertion.eq(set(NAMES), qmflows._REQUIRES_RDKIT.keys())


@pytest.mark.skipif(sys.version_info < (3, 7), reason="requires python 3.7")
@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_namespace() -> None:
    """Test that ``qmflows._REQUIRES_RDKIT`` and the sub-modules' ``__all__`` are synced."""
    name_set = {"components", "examples"}
    name_set.update(qmflows.components.__all__)
    name_set.update(qmflows.examples.__all__)

    assertion.eq(name_set, qmflows._REQUIRES_RDKIT.keys())
