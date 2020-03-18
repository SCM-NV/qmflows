"""Interface to call the ADF KFReader."""
__all__ = ['kfreader']

from scm import plams


def kfreader(path_t21, section=None, prop=None):
    """Use the plams KFfile to read the TAPE21 File."""
    kf = plams.tools.kftools.KFFile(path_t21)
    return kf.read(section, prop)
