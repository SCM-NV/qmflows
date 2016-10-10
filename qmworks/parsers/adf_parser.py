
import plams


def kfreader(path_t21, section=None, prop=None):
    """
    Use the plams KFfile to read the TAPE21 File.
    """
    kf = plams.kftools.KFFile(path_t21)
    return kf.read(section, prop)
