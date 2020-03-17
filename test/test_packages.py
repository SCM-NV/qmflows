"""Testing the workflow runner."""
from assertionlib import assertion
from noodles import schedule

from qmflows import run


@schedule
def f(a, b):
    return a * b


@schedule
def g(c, d):
    return 2 * c + d


def test_run_packages(tmpdir):
    """Test the Workflow runner."""
    wf = g(3, f(5, 2))
    result = run(wf, path=tmpdir, folder=tmpdir)
    assertion.eq(result, 16)
