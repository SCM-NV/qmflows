
from qmworks.settings import Settings
import plams
import numpy as np


def kfreader(path_t21, section=None, prop=None):
    """
    Use the plams KFfile to read the TAPE21 File.
    """
    kf = plams.tools.kftools.KFFile(path_t21)
    return kf.read(section, prop)


def extract_properties_rkf(path_rkf, key=None):
    """
    Read result from a DFTB computation using the job_name.rkf file.
    """
    kf = plams.tools.kftools.KFFile(path_rkf).read
    props = Settings()

    for i in range(kf('Properties', 'nEntries')):
        typ = kf('Properties', 'Type(' + str(i + 1) + ')').strip()
        subtype = kf('Properties', 'Subtype(' + str(i + 1) + ')').strip()
        value = kf('Properties', 'Value(' + str(i + 1) + ')')
        props[typ][subtype] = value
    ret = props[key]
    if isinstance(ret, list):
        ret = np.array(ret)
    return ret
