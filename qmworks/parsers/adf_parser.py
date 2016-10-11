
from qmworks import Settings
import plams


def kfreader(path_t21, section=None, prop=None):
    """
    Use the plams KFfile to read the TAPE21 File.
    """
    kf = plams.kftools.KFFile(path_t21)
    return kf.read(section, prop)


def extract_properties_rkf(path_rkf, key=None):
    """
    """
    kf = plams.kftools.KFFile(path_rkf)
    props = Settings()
    
    for i in range(kf('Properties', 'nEntries')):
        typ = kf('Properties', 'Type(' + str(i + 1) + ')').strip()
        subtype = kf('Properties', 'Subtype(' + str(i + 1) + ')').strip()
        value = kf('Properties', 'Value(' + str(i + 1) + ')')
        props[typ][subtype] = value

    return props[key]
