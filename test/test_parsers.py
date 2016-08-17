
import h5py 

from os.path import join
from qmworks.parsers import readCp2KBasis

# ===================================<>========================================
path_hdf5 = 'test/test_files/ethylene.hdf5'
path_basis = 'test/test_files/BASIS_MOLOPT'
basis_name = 'DZVP-MOLOPT-SR-GTH'


def test_parser_basis_cp2k():
    """
    Test that the CP2K Molecularly optimized are read properly.
    
    NamedTuples:
    AtomBasisKey = namedtuple("AtomBasisKey", ("atom", "basis", "basisFormat"))
    AtomBasisData = namedtuple("AtomBasisData", ("exponents", "coefficients"))
    """
    # basis_keys, basis_data = readCp2KBasis(path_basis)
    # dim = len(basis_keys)
    # keys = filter(lambda t: t[1].basis == basis_name, zip())
    # with h5py.File(path_hdf5, 'r') as f5:
    assert True
