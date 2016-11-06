
__all__ = ['StoreasHDF5', 'cp2k2hdf5', 'read_from_hdf5',
           'turbomole2hdf5']

# ==========> Standard libraries and third-party <===============
from functools import partial
from os.path import join

import h5py
import numpy as np
# ==================> Internal modules <==========
from qmworks.parsers.cp2KParser import readCp2KBasis
from qmworks.parsers.turbomoleParser import readTurbomoleBasis

# ====================><==============================


def read_from_hdf5(path_hdf5, path_to_properties):
    try:
        with h5py.File(path_hdf5, 'r') as f5:
            if isinstance(path_to_properties, list):
                return [f5[path].value for path in path_to_properties]
            else:
                return f5[path_to_properties].value
    except KeyError:
        msg = "There is not {} stored in the HDF5\n".format(path_to_properties)
        raise KeyError(msg)
    except FileNotFoundError:
        msg = "there is not HDF5 file containing the numerical results"
        raise RuntimeError(msg)


class StoreasHDF5:
    """
    Class to store inside a HDF5 file numerical array with optional attributes.
    """
    def __init__(self, file_h5, packageName):
        self.file_h5 = file_h5
        self.name = packageName

    def funHDF5(self, pathProperty, data):
        """
        creates a data set using ``data`` and saves the data using
        ``pathProperty`` in the HDF5 file.

        :param pathProperty: path to store the property in HDF5.
        :type pathProperty: String
        :param data: Numeric array containing the property.
        :type data: Numpy array
        :returns: **None**
        """
        self.file_h5.require_dataset(pathProperty, shape=np.shape(data),
                                     data=data, dtype=np.float32)

    def funHDF5_attrs(self, nameAttr, attr, pathProperty, data):
        """
        creates a data set using ``data`` and some attributes.

        :param nameAttr: Name of the attribute assoaciated with the data.
        :type nameAttr: String
        :param attr: Actual atttribute.
        :type attr: String | Numpy array
        :param pathProperty: path to store the property in HDF5.
        :type pathProperty: String
        :param data: Numeric array containing the property.
        :type data: Numpy array
        :returns: **None**
        """

        dset = self.file_h5.require_dataset(pathProperty, shape=np.shape(data),
                                            data=data, dtype=np.float32)
        dset.attrs[nameAttr] = attr

    def saveBasis(self, parserFun, pathBasis):
        """
        Store the basis set.

        :param parserFun: Function to parse the file containing the
                          information about the primitive contracted Gauss
                          functions.
        :param pathBasis: Absolute path to the file containing the basis
                          sets information.
        :type pathBasis: String.
        :returns: **None**
        """

        keys, vals = parserFun(pathBasis)
        pathsExpo = [join(self.name, "basis", xs.atom, xs.basis, "exponents")
                     for xs in keys]
        pathsCoeff = [join(self.name, "basis", xs.atom, xs.basis,
                           "coefficients") for xs in keys]

        for ps, es in zip(pathsExpo, [xs.exponents for xs in vals]):
            self.funHDF5(ps, es)

        fss = [xs.basisFormat for xs in keys]
        css = [xs.coefficients for xs in vals]

        # save basis set coefficients and their correspoding format
        for path, fs, css in zip(pathsCoeff, fss, css):
            self.funHDF5_attrs("basisFormat", str(fs), path, css)


# =======================================================
# CP2K Interface

def cp2k2hdf5(file_h5, keys):
    """
    Use a list of namedtuple ``keys`` to retrieve information
    from text otuput files and store it in HDF5 format
    """
    for k in keys:
        cp2kOpts(file_h5, k)


def cp2kOpts(file_h5, key):
    """
    Read from a text file some numerical information and store it in HDF5
    Format. The available options to read and store are:
    - basis set
    - Molecular orbitals
    -Overlap Mtrix
    """
    storeCp2k = StoreasHDF5(file_h5, "cp2k")

    args = key.args
    name = key.name
    d = {"basis": partial(storeCp2k.saveBasis, readCp2KBasis)}

    return d[name](*args)


# ==================================================
# TurboMol Interface


def turbomole2hdf5(file_h5, keys):
    for k in keys:
        turbomoleOpts(file_h5, k)


def turbomoleOpts(file_h5, key):
    """
    Read from text files numerical properties and store them in HDF5 format.
    Available options:
    - Basis set
    - Molecular Orbitals
    """
    storeTurbo = StoreasHDF5(file_h5, "turbomole")
    args = key.args
    name = key.name
    d = {"basis": partial(storeTurbo.saveBasis, readTurbomoleBasis)}

    return d[name](*args)
