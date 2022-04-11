"""A pytest ``conftest.py`` file."""

import os
from pathlib import Path

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
