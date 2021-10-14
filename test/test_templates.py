"""Test that the templates behave correctly."""

import sys
import qmflows
from qmflows import Settings, templates
from qmflows.templates import freq, geometry, singlepoint, ts
from assertionlib import assertion

import pytest


def transform(x: Settings) -> Settings:
    """Transform `x` to dict and back to Settings."""
    return Settings(x.as_dict())


def test_templates():
    """Test that the yaml files are read properly."""
    b1 = freq == transform(freq)
    b2 = geometry == transform(geometry)
    b3 = singlepoint == transform(singlepoint)
    b4 = ts == transform(ts)

    assertion.truth(all((b1, b2, b3, b4)))


@pytest.mark.parametrize("name", templates.__all__)
@pytest.mark.skipif(sys.version_info < (3, 7), reason="Requires python >= 3.7")
def test_id(name: str) -> None:
    """Test that getting a template returns a copy."""
    s1 = getattr(qmflows, name)
    s2 = getattr(qmflows, name)
    assertion.eq(s1, s2)
    assertion.is_not(s1, s2)


def test_namespace() -> None:
    assertion.issubset(templates.__all__, dir(qmflows))
