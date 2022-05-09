"""Test the :mod:`sphinx` documentation generation."""

import warnings
from pathlib import Path

import pytest
from qmflows.test_utils import stdout_to_logger

try:
    from sphinx.application import Sphinx
    from sphinx.errors import SphinxWarning
except ImportError:
    HAS_SPHINX = False
else:
    HAS_SPHINX = True

SRCDIR = CONFDIR = 'docs'


@pytest.mark.skipif(not HAS_SPHINX, reason='Requires Sphinx')
@pytest.mark.xfail
def test_sphinx_build(tmp_path: Path) -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    try:
        app = Sphinx(
            SRCDIR,
            CONFDIR,
            tmp_path / "build",
            tmp_path / "build" / "doctrees",
            buildername='html',
            warningiserror=True,
            status=stdout_to_logger,
        )
        app.build(force_all=True)
    except SphinxWarning as ex:
        # Do not raise if the exception is connection related
        if "Max retries exceeded" not in str(ex):
            raise
        else:
            warning = RuntimeWarning(str(ex))
            warning.__cause__ = ex
            warnings.warn(warning)
            pytest.xfail(str(ex))
