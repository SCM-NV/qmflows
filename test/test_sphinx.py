"""Test the :mod:`sphinx` documentation generation."""

import shutil
import warnings
from os.path import join, isdir

from sphinx.application import Sphinx
from sphinx.errors import SphinxWarning

SRCDIR = CONFDIR = 'docs'
OUTDIR = join('tests', 'test_files', 'build')
DOCTREEDIR = join('tests', 'test_files', 'build', 'doctrees')


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
            warnings.warn(ex)
    finally:
        if isdir(OUTDIR):
            shutil.rmtree(OUTDIR)
