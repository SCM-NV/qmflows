

__all__ = ['lookup_H5']

# ==========> Standard libraries and third-party <===============
from os.path import (dirname, join)

# ====================<>==============================
# HDF5 Functions


def lookup_H5(key, file_H5):

    def lookup(name, obj):
        dir_path = dirname(name)
        key_path = join(dir_path, key)

        if name == key_path:
            return obj
        else:
            return None

    return file_H5.visititems(lookup)


def existGroup(file_H5, key):

    def lookup(name, obj):
        dir_path = dirname(name)
        key_path = join(dir_path, key)

        if name == key_path:
            return True
        else:
            return False

    return file_H5.visititems(lookup)
