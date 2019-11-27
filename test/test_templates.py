from qmflows.templates import (freq, geometry, singlepoint, ts)
from qmflows.utils import (dict2Setting, settings2Dict)


def fun(x):
    """Assert that the functions are the inverse of each other."""
    return dict2Setting(settings2Dict(x))


def test_templates():
    """
    Test that the JSON files are read properly.
    """
    b1 = freq == fun(freq)
    b2 = geometry == fun(geometry)
    b3 = singlepoint == fun(singlepoint)
    b4 = ts == fun(ts)

    assert all([b1, b2, b3, b4])
