
from os.path import join
from qmworks.common import InputKey
from qmworks.hdf5.quantumHDF5 import turbomole2hdf5

import h5py
import os
# ====================================<>=======================================
path_hdf5 = 'test/test_files/test.hdf5'


def test_store_basisSet():
    """
    Check if the turbomole basis set are read
    and store in HDF5 format.
    """
    path_basis = 'test/test_files/basis_turbomole'
    keyBasis = InputKey("basis", [path_basis])

    try_to_remove(path_hdf5)
    with h5py.File(path_hdf5, chunks=True) as f5:
        try:
            # Store the basis sets
            turbomole2hdf5(f5, [keyBasis])
            if not f5["turbomole/basis"]:
                assert False
        except RuntimeError:
            try_to_remove(path_hdf5)
            assert False

            
def test_store_MO_h5():
    """
    test if the MO are stored in the HDF5 format
    """
    path = join('/turbomole', 'test', 'ethylene')
    path_es = join(path, 'eigenvalues')
    path_css = join(path, 'coefficients')
    number_of_orbs = 36
    number_of_orb_funs = 38

    try_to_remove(path_hdf5)
    with h5py.File(path_hdf5) as f5:
        path_es, path_css = dump_MOs_coeff(f5, path_es, path_css,
                                           number_of_orbs, number_of_orb_funs)
        if f5[path_es] and f5[path_css]:
            try_to_remove(path_hdf5)
            assert True
        else:
            try_to_remove(path_hdf5)
            assert False
            
            
def dump_MOs_coeff(handle_hdf5, path_es, path_css, number_of_orbs,
                   number_of_orb_funs):
    """
    MO coefficients are stored in row-major order, they must be transposed
    to get the standard MO matrix.
    :param files: Files to calculate the MO coefficients
    :type  files: Namedtuple (fileXYZ,fileInput,fileOutput)
    :param job: Output File
    :type  job: String
    """
    key = InputKey('orbitals', [path_MO, number_of_orbs, number_of_orb_funs,
                                path_es, path_css])

    turbomole2hdf5(handle_hdf5, [key])

    return path_es, path_css


            
def try_to_remove(path):
    """
    Remove a file if it exists
    """
    try:
        os.remove(path)
    except OSError:
        pass
            
