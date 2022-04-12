import pytest
from assertionlib import assertion
from packaging.version import Version

from qmflows import version_info, __version__


def test_pep_440() -> None:
    assertion.assert_(Version, __version__)


@pytest.mark.parametrize("name", ["major", "minor", "micro", "patch", "bug", "maintenance"])
def test_version_info(name: str) -> None:
    attr: int = getattr(version_info, name)
    assertion.isinstance(attr, int)
    assertion.ge(attr, 0)
