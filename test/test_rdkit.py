"""Tests related to the optional rdkit dependency."""

import sys
import importlib

import pytest
import qmflows
from qmflows.test_utils import HAS_RDKIT
from qmflows._init_utils import RDKIT_SET
from assertionlib import assertion


@pytest.mark.parametrize("name", sorted(RDKIT_SET))
def test_sub_module(name: str) -> None:
    """Test :func:`getattr` operations on rdkit-requiring objects."""
    if HAS_RDKIT:
        assert getattr(qmflows, name)
    else:
        match = f"{name!r} requires the optional RDKit package"
        with pytest.raises(ImportError, match=match):
            getattr(qmflows, name)


@pytest.mark.skipif(sys.version_info < (3, 7), reason="requires python 3.7")
@pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
def test_rdkit() -> None:
    """Test that ``qmflows._init_utils.RDKIT_DICT`` and the sub-modules' ``__all__`` are synced."""
    name_set = {"components", "examples"}
    name_set.update(qmflows.components.__all__)
    name_set.update(qmflows.examples.__all__)
    assertion.eq(name_set, RDKIT_SET)


def test_dir() -> None:
    """Test that RDKit functions are in-/excluded from ``dir``."""
    all_names = RDKIT_SET - {"components", "examples"}
    if HAS_RDKIT:
        assertion.issubset(RDKIT_SET, dir(qmflows))
        assertion.issubset(all_names, qmflows.__all__)
    else:
        assertion.isdisjoint(RDKIT_SET, dir(qmflows))
        assertion.isdisjoint(all_names, qmflows.__all__)
