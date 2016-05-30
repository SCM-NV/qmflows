
__author__ = "Felipe Zapata"

__all__ = ['loadFromJSON', 'json2Settings', 'filter_path',
           'or_pattern', 'find_tape21', 'search_environ_var',
           'find_cp2k_basis', 'find_cp2k_MO']

# ================> Python Standard  and third-party <==========
import fnmatch
import json
import os
import pkg_resources as pkg

from os.path import join

# ==================> Internal modules <==========
from qmworks.utils import dict2Setting
# ======================> <=========================

# ================> JSON Utilities <==========================


def loadFromJSON(fileName):
    """
    """
    with open(fileName, 'r') as f:
        xs = json.load(f)
    return xs


# def json2Settings(file_JSON):
#     xs = loadFromJSON(file_JSON)
#     return dict2Setting(xs)


def json2Settings(xs):
    """
    transform a string containing some data in JSON format
    to a Settings object
    """
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = json.loads(xs)  # Json object must be string
    return dict2Setting(s)


# ==================> Files Utilities <==================


def filter_path(path, pattern):
    files = os.listdir(path)
    ps = fnmatch.filter(files, pattern)
    return [os.path.join(path, x) for x in ps]


def or_pattern(path, pattern1, pattern2):
    files = os.listdir(path)
    xs = fnmatch.filter(files, pattern1)
    if not xs:
        xs = fnmatch(files, pattern2)
    return [os.path.join(path, x) for x in xs]


def find_tape21(path):
    xs = or_pattern(path, "*.t21", "*TAPE21*")
    if not xs:
        err = 'There is not a tape21 file in path:{} '.format(path)
        raise FileError(err)
    else:
        return xs


def search_environ_var(var):
    """
    Looks if the environmental variable ``var`` is defined.
    """
    try:
        defaultPath = os.environ[var]
    except KeyError:
        pass

    if not defaultPath:
        msg = 'There is not an environmental variable called: {}'.format(var)
        raise EnvironmentError(msg)

    return defaultPath


def find_cp2k_basis(path):
    """
    """
    try:
        path_basis = os.environ['BASISCP2K']
    except KeyError:
        path_basis = filter_path(path, "*[Bb][Aa][Ss][Ii][Ss]*")

    if not path_basis:
        msg = 'I have not found a file containing the Cp2Kbasis.'
        'You can provide the path to this file using the following unix'
        'enviromental variable: BASISCP2K=Path/To/The/CP2K/Basis/Set'
        raise FileNotFoundError(msg)
    return path_basis


def find_cp2k_MO(path):
    xs = or_pattern(path, "mo*Log", "mo*out")
    print("Mo_files: ")
    if not xs:
        msg = 'There is not an output file containing the cp2k MO in path:{}\n'
        err = msg.format(path)
        raise FileNotFoundError(err)
    else:
        return xs[0]
