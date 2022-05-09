"""A pytest ``conftest.py`` file."""

import os
import sys
import types
import tempfile
import importlib
import contextlib
from typing import Generator
from pathlib import Path

import pytest
from scm.plams import config

from qmflows import InitRestart

_ROOT = Path("src") / "qmflows"
_collect_ignore = [
    _ROOT/ "settings.py",
    _ROOT / "parsers" / "xyzParser.py",
    _ROOT / "parsers" / "generic_parsers.py",
    _ROOT / "parsers" / "parser.py",
    _ROOT / "parsers" / "orca_parser.py",
    _ROOT / "parsers" / "cp2KParser.py",
    _ROOT / "parsers" / "adf_parser.py",
    _ROOT / "templates" / "templates.py",
    _ROOT / "packages" / "cp2k_mm.py",
    _ROOT / "packages" / "cp2k_package.py",
    _ROOT / "packages" / "orca.py",
    _ROOT / "packages" / "package_wrapper.py",
    _ROOT / "packages" / "packages.py",
    _ROOT / "packages" / "serializer.py",
    _ROOT / "packages" / "SCM.py",
]

collect_ignore = [str(i) for i in _collect_ignore]


def _del_all_attr(module: types.ModuleType) -> None:
    """Delete all module attributes so they will be reloaded upon reloading the module."""
    attr_names = [i for i in vars(module) if not (i.startswith("__") and i.endswith("__"))]
    for name in attr_names:
        delattr(module, name)


@pytest.fixture(autouse=True, scope="session")
def reload_qmflows() -> Generator[None, None, None]:
    """Reload qmflows in order to re-trigger coverage.

    This is necasary as the setup.cfg ``filterwarnings`` option will load qmflows prior
    to code coverage being turned on, thus falsely excluding a large portion of qmflows
    from the code coverage report.
    s
    """
    yield None

    module_names = [i for i in sys.modules if i.startswith("qmflows")]
    for name in module_names:
        module = sys.modules.pop(name)
        _del_all_attr(module)

    importlib.import_module("qmflows")


@pytest.fixture(autouse=True, scope="session")
def configure_plams_logger() -> "Generator[None, None, None]":
    """Remove the date/time prefix from the PLAMS logging output."""
    # Ensure the plams.config dict is populated by firing up plams.init once
    with open(os.devnull, "w") as f1, tempfile.TemporaryDirectory() as f2:
        with contextlib.redirect_stdout(f1), InitRestart(f2):
            pass

    assert "log" in config
    log_backup = config.log.copy()
    config.log.time = False
    config.log.date = False

    yield None
    config.log = log_backup
