from noodles import schedule
from qmflows import run
import os
import shutil
import tempfile


@schedule
def f(a, b):
    return a * b


@schedule
def g(c, d):
    return 2 * c + d


def test_run_packages():
    """
    Test the Workflow runner.
    """
    folder = tempfile.mkdtemp(prefix="qmflows_")
    try:
        wf = g(3, f(5, 2))
        result = run(wf, path=folder, folder=folder)
        assert result == 16
    finally:
        if os.path.isdir(folder):
            shutil.rmtree(folder)
