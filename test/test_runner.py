"""Testing the workflow runner."""

from pathlib import Path

from assertionlib import assertion
from noodles import schedule

from qmflows import run


@schedule
def f(a: int, b: int) -> int:
    return a * b


@schedule
def g(c: int, d: int) -> int:
    return 2 * c + d


def test_run_packages(tmp_path: Path) -> None:
    """Test the Workflow runner."""
    wf = g(3, f(5, 2))
    result = run(wf, path=tmp_path, folder="test_run_packages")
    assertion.eq(result, 16)
