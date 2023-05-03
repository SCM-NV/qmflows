"""Interface to call the ADF KFReader."""

from __future__ import annotations

from typing import Any

from scm import plams

from ..type_hints import PathLike

__all__ = ['kfreader']


def kfreader(path_t21: PathLike, section: None | str = None,
             prop: None | str = None) -> Any:
    """Use the plams KFfile to read the TAPE21 File."""
    kf = plams.tools.kftools.KFFile(path_t21)
    return kf.read(section, prop)
