from noodles import (get_workflow, schedule_hint)
from qmflows.draw_workflow import draw_workflow
import os
import pytest
import tempfile


@schedule_hint()
def f(a, b):
    return a ** b


@schedule_hint()
def g(c, d):
    return c + d * 2


@pytest.mark.xfail
def test_draw_workflow():
    filename = tempfile.mktemp(suffix='.svg', prefix='workflow_')

    try:
        wf = f(3, g(5, 2))
        draw_workflow(filename, get_workflow(wf))
        assert os.path.isfile(filename)
    except:
        assert False
    finally:
        if os.path.isfile(filename):
            os.remove(filename)
