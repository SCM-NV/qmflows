from nose.plugins.attrib import attr
from qmflows.common import InputKey
from qmflows.hdf5.quantumHDF5 import turbomole2hdf5

import h5py
import os
# ====================================<>=======================================
path_hdf5 = 'test/test_files/test.hdf5'


@attr('fast')
def test_store_basisSet():
    """
    Check if the turbomole basis set are read
    and store in HDF5 format.
    """
    path_basis = 'test/test_files/basis_turbomole'
    keyBasis = InputKey("basis", [path_basis])

    with h5py.File(path_hdf5) as f5:
        try:
            # Store the basis sets
            turbomole2hdf5(f5, [keyBasis])
            if not f5["turbomole/basis"]:
                assert False
        finally:
            try_to_remove(path_hdf5)


def try_to_remove(path):
    """
    Remove a file if it exists
    """
    try:
        os.remove(path)
    except OSError:
        pass
