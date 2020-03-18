"""Test that the templates behave correctly."""
from qmflows import Settings
from qmflows.templates import freq, geometry, singlepoint, ts
from assertionlib import assertion


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
