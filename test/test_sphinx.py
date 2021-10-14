"""Test the :mod:`sphinx` documentation generation."""

import shutil
import warnings
from os.path import join, isdir

import pytest

try:
    from sphinx.application import Sphinx
    from sphinx.errors import SphinxWarning
except ImportError:
    HAS_SPHINX = False
else:
    HAS_SPHINX = True

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('test', 'test_files', 'build')
DOCTREEDIR = join('test', 'test_files', 'build', 'doctrees')


@pytest.mark.skipif(not HAS_SPHINX, reason='Requires Sphinx')
def test_sphinx_build() -> None:
    """Test :meth:`sphinx.application.Sphinx.build`."""
    try:
        app = Sphinx(SRCDIR, CONFDIR, OUTDIR, DOCTREEDIR,
                     buildername='html', warningiserror=True)
        app.build(force_all=True)
    except SphinxWarning as ex:
        # Do not raise if the exception is connection related
        if "Max retries exceeded" not in str(ex):
            raise
        else:
            warning = RuntimeWarning(str(ex))
            warning.__cause__ = ex
            warnings.warn(warning)
    finally:
        if isdir(OUTDIR):
            shutil.rmtree(OUTDIR)
