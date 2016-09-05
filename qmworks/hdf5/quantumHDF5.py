
__author__ = "Felipe Zapata"

__all__ = ['StoreasHDF5', 'adf2hdf5', 'save_ADF_Info', 'cp2k2hdf5', 'turbomole2hdf5']


# ==========> Standard libraries and third-party <===============
import numpy as np
from collections import namedtuple
from functools import partial
from os.path import join

# ==================> Internal modules <==========
from plams import KFReader
from qmworks.common import InfoMO
from qmworks.parsers.cp2KParser import (readCp2KBasis, readCp2KCoeff, readCp2KOverlap)
from qmworks.parsers.turbomoleParser import (readTurbomoleBasis, readTurbomoleMO)
from qmworks.utils import zipWith


# ====================><==============================

class StoreasHDF5:
    """
    Class to store inside a HDF5 file numerical array with optional attributes.
    """
    def  __init__(self, file_h5, packageName):
        self.file_h5 = file_h5
        self.name = packageName

    def funHDF5(self, pathProperty, data):
        """
        creates a data set using ``data`` and saves the data using ``pathProperty``
        in the HDF5 file.

        :parameter pathProperty: path to store the property in HDF5.
        :type pathProperty: String
        :parameter data: Numeric array containing the property.
        :type data: Numpy array
        :returns: **None**
        """
        self.file_h5.require_dataset(pathProperty, shape=np.shape(data), data=data,
                                     dtype=np.float32)

    def funHDF5_attrs(self, nameAttr, attr, pathProperty, data):
        """
        creates a data set using ``data`` and some attributes.

        :parameter nameAttr: Name of the attribute assoaciated with the data.
        :type nameAttr: String
        :parameter attr: Actual atttribute.
        :type attr: String | Numpy array
        :parameter pathProperty: path to store the property in HDF5.
        :type pathProperty: String
        :parameter data: Numeric array containing the property.
        :type data: Numpy array
        :returns: **None**
        """

        dset = self.file_h5.require_dataset(pathProperty, shape=np.shape(data),
                                            data=data, dtype=np.float32)
        dset.attrs[nameAttr] = attr

    def saveBasis(self, parserFun, pathBasis):
        """
        Store the basis set.

        :parameter parserFun: Function to parse the file containing the
        information about the primitive contracted Gauss functions.
        :parameter pathBasis: Absolute path to the file containing the basis
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

    def saveMO(self, parserFun, pathMO, nOrbitals=None, nOrbFuns=None,
               pathEs=None, pathCs=None, nOccupied=None, nHOMOS=None,
               nLUMOS=None):
        """
        Save Molecular orbital eigenvalues and eigenvectors.

        :parameter parserFun: Function to parse the file containing the MOs.
        :parameter pathMO: Absolute path to the MOs file.
        :type pathMO: String
        :parameter nOrbitals: Number of orbitals
        :type nOrbitals: Int
        :param nOrbFuns: Number of Atomic orbitals function that made up
        each of the molecular orbitals.
        :type nOrbFuns: Int
        :parameter pathEs: Path to the MO eigenvalues in the HDF5 file
                           (default is: <softwareName>/mo/eigenvalues).
        :type pathEs: String
        :parameter pathCs: Path to the MO coefficients in the HDF5 file
                           (default is: <softwareName>/mo/coefficients).
        :type pathCs: String
        :param nOccupied: Number of occupied molecular orbitals.
        :type nOccupied: Int
        :param nHOMOS: number of HOMO orbitals to store in HDF5.
        :type nHOMOS: Int
        :param nLUMOS: number of HUMO orbitals to store in HDF5.
        :type nLUMOS: Int
        :returns: **None**
        """

        pathEs = pathEs if pathEs is not None else join(self.name, "mo",
                                                        "eigenvalues")
        pathCs = pathCs if pathCs is not None else join(self.name, "mo",
                                                        "coefficients")

        # Save all the frontier orbitals.
        # NOTE: IT ASSUMES THAT The USER HAVE SELECTED A RANGE OF MO
        # TO PRINT.
        if not (nHOMOS is None and nLUMOS is None):
            infoMO = parserFun(pathMO, nHOMOS + nLUMOS, nOrbFuns)

        elif nOrbitals is not None:
            infoMO = parserFun(pathMO, nOrbitals, nOrbFuns)
            if nOrbitals > nHOMOS + nLUMOS:
                # Drop Coefficients that below and above nHOMOS and nLUMOS, respectively.
                ess, css  = infoMO
                css = np.transpose(css)
                eigenVals = ess[nOccupied - nHOMOS: nOccupied + nLUMOS]
                coefficients = css[nOccupied - nHOMOS: nOccupied + nLUMOS]
                infoMO = InfoMO(eigenVals, np.transpose(coefficients))

        zipWith(self.funHDF5)([pathEs, pathCs])([infoMO.eigenVals, infoMO.coeffs])

    def saveOverlap(self, parserFun, pathOverlap, nOrbitals=None, path=None):
        """
        Store the Overlap matrix.

        :parameter parserFun: Function to parse the file containing the Overlap
                              matrix.
        :parameter pathOverlap: Absolute path to the overlap matrix file.
        :type pathOverlap: String
        :parameter nOrbitals: Number of orbitals
        :type nOrbitals: Int
        :parameter path: Path to store the Overlap matrix in the HDF5 file.
        :type path: String
        :returns: **None**
        """
        mtx = parserFun(pathOverlap, nOrbitals)
        path = path if path is not None  else join(self.name, "overlap")
        self.funHDF5(path, mtx)

# ====================><==============================
# named Tuples

H5Args = namedtuple("H5Args", ("path", "h5File", "keys"))


# Quantum package output dumping to HDF5
# ADF Interfaces

def adf2hdf5(file_h5, keys):
    """
    serialize the numerical arrays specfied by ``keys`` named tuples.

    :parameter file_h5: HDF5 file storing the numerical arrays
    :type      file_h5: h5py handler
    :parameter    keys: Properties to store
    :type         keys: Namedtuple
    """
    for k in keys:
        adfOpts(file_h5, k)


def adfOpts(file_h5, key):
    """
    Store the ``key`` properties in the ``file_h5``
    """
    path = key.args
    (section, value) = dict_genericKey2ADF[key.name]
    kfreader = KFReader(path)
    save_ADF_Info(kfreader, section, value, file_h5)


def save_ADF_Info(kfreader, section, value, file_h5):
    """
    Saves in the HDF5 file ``value`` from the ``section``.
    """
    path = join(section, value)
    xs = kfreader.read(section, value)
    arr = np.array(xs)

    file_h5.create_dataset(path, data=arr, dtype=np.float32)

dict_genericKey2ADF = {'energy': (), 'orbitals': (),
                       'gradient': ('Geopt', 'Gradients'),
                       'hessian': ('Freq', 'Hessian_CART')}
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
    d = {"basis": partial(storeCp2k.saveBasis, readCp2KBasis),
         "orbitals": partial(storeCp2k.saveMO, readCp2KCoeff),
         "overlap": partial(storeCp2k.saveOverlap, readCp2KOverlap)}

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
    d = {"basis": partial(storeTurbo.saveBasis, readTurbomoleBasis),
         "orbitals": partial(storeTurbo.saveMO, readTurbomoleMO)}

    return d[name](*args)

# =========================================================
# Dirac Interface


def dirac2hdf5():
    pass

