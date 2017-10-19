from nose.plugins.attrib import attr
from qmflows.templates import (freq, geometry, singlepoint, ts)
from qmflows.utils import (dict2Setting, settings2Dict)


@attr('fast')
def test_templates():
    """
    Test that the JSON files are read properly.
    """
    fun = lambda x: dict2Setting(settings2Dict(x))
    b1 = freq == fun(freq)
    b2 = geometry == fun(geometry)
    b3 = singlepoint == fun(singlepoint)
    b4 = ts == fun(ts)

    assert all([b1, b2, b3, b4])
