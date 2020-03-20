"""Interface to call the ADF KFReader."""
__all__ = ['kfreader']

from typing import Optional, Any

from scm import plams

from ..type_hints import PathLike


def kfreader(path_t21: PathLike, section: Optional[str] = None,
             prop: Optional[str] = None) -> Any:
    """Use the plams KFfile to read the TAPE21 File."""
    kf = plams.tools.kftools.KFFile(path_t21)
    return kf.read(section, prop)
