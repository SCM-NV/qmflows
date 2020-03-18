"""Test for PEP8 compliance."""

from typing import Tuple
import os
import textwrap

import pycodestyle  # formerly known as pep8

INCLUDE_PATHS: Tuple[str, ...] = ('src',)
EXCLUDE_PATHS: Tuple[str, ...] = ()


def test_pep8_conformance() -> None:
    """Test for PEP8 compliance."""
    print(f"PEP8 check of directories: {', '.join(INCLUDE_PATHS)}\n")

    # Get paths wrt package root
    root = os.path.dirname(os.path.dirname(__file__))
    include = [os.path.join(root, path) for path in INCLUDE_PATHS]
    exclude = [os.path.join(root, path) for path in EXCLUDE_PATHS]

    # Increase the maximum amount of characters per line from 79 to 100
    style = pycodestyle.StyleGuide(max_line_length=100)
    style.options.exclude.extend(exclude)

    success = style.check_files(include).total_errors == 0

    if not success:
        print(textwrap.dedent("""
            Your Python code does not conform to the official Python style
            guide (PEP8), see https://www.python.org/dev/peps/pep-0008

            A list of warning and error messages can be found above,
            prefixed with filename:line number:column number.

            Run `yapf -i yourfile.py` to automatically fix most errors.
            Run `yapf -d yourfile.py` to preview what would be changed.
            Run `pip install --upgrade yapf` to install the latest version
            of yapf.
        """))

    assert success, "Your code does not conform to PEP8"
