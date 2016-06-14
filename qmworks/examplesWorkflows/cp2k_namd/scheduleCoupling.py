__author__ = "Felipe Zapata"

# ================> Python Standard  and third-party <==========
from noodles import gather, schedule  # Workflow Engine
from os.path import join

import h5py
import os
import numpy as np
# ==================> Internal modules <==========
from nac.integrals import calc_transf_matrix, calculateCoupling3Points

from qmworks.common import AtomXYZ
from qmworks.hdf5.quantumHDF5 import StoreasHDF5
from qmworks.parsers import parse_string_xyz

# ================> Constants <================
angs2au = 1 / 0.529177249  # Angstrom to a.u
femtosec2au = 1 / 2.41888432e-2  # from femtoseconds to au

# ==============================> Schedule Tasks <=========================


@schedule
def schedule_transf_matrix(path_hdf5, atoms, basisName, work_dir, packageName):
    """
    calculate the transformation of the overlap matrix from both spherical
    to cartesian and from cartesian to spherical.

    :returns: Numpy matrix containing the transformation matrix.
    """
    with h5py.File(path_hdf5) as f5:
        mtx = calc_transf_matrix(f5, atoms, basisName, packageName)
        store = StoreasHDF5(f5, packageName)
        path = os.path.join(work_dir, 'trans_mtx')
        store.funHDF5(path, mtx)
    return path


@schedule
def lazy_schedule_couplings(i, path_hdf5, dictCGFs, geometries, mo_paths, dt=1,
                            hdf5_trans_mtx=None, units='angstrom',
                            output_folder=None):
    """
    :param i: nth coupling calculation
    :type i: int
    :paramter dictCGFS: Dictionary from Atomic Label to basis set
    :type     dictCGFS: Dict String [CGF],
              CGF = ([Primitives], AngularMomentum),
              Primitive = (Coefficient, Exponent)
    :parameter geometries: Tuple molecular geometries stored as strings
    :type      geometries: str, str, str)
    :parameter coefficients: Tuple of Molecular Orbital coefficients.
    :type      coefficients: (Matrix, Matrix, Matrix)
    :parameter dt: dynamic integration time
    :type      dt: Float (Femtoseconds)
    """
    def calc_coupling(output_path, dt):
    
        if hdf5_trans_mtx is not None:
            trans_mtx = retrieve_hdf5_data(path_hdf5, hdf5_trans_mtx)
        else:
            trans_mtx = None

        xss = tuple(map(parse_string_xyz, geometries))

        if 'angstrom' in units.lower():
            xss = tuple(map(change_mol_units, xss))
            
        dt_au = dt * femtosec2au
        
        mos = tuple(map(lambda j: retrieve_hdf5_data(path_hdf5, mo_paths[i + j][1]), range(3)))

        rs = calculateCoupling3Points(xss, mos, dictCGFs, dt_au, trans_mtx)
        
        with h5py.File(path_hdf5) as f5:
            store = StoreasHDF5(f5, 'cp2k')
            store.funHDF5(output_path, rs)

    if output_folder is None:
        msg = 'There was not specified a path in the HDF5 file to store the coupling\n'
        raise RuntimeError(msg)

    output_path = join(output_folder, 'coupling_{}'.format(i))

    try:
        with h5py.File(path_hdf5, 'r') as f5:
            f5[output_path]
    except KeyError:
        # If the coupling is not store in the HDF5 calculate it
        calc_coupling(output_path, dt)

    return output_path


def write_hamiltonians(path_hdf5, work_dir, mo_paths, path_couplings,
                       nPoints, path_dir_results=None):
    """
    Write the real and imaginary components of the hamiltonian using both
    the orbitals energies and the derivative coupling accoring to:
    http://pubs.acs.org/doi/abs/10.1021/ct400641n
    **Units are: Rydbergs**.
    """
    def write_pyxaid_format(arr, fileName):
        np.savetxt(fileName, arr, fmt='%10.5e', delimiter='  ')

    ham_files = []
    for i in range(nPoints):
        path_coupling = path_couplings[i]
        css = retrieve_hdf5_data(path_hdf5, path_coupling)
        energies = retrieve_hdf5_data(path_hdf5, mo_paths[i][0])
            
        file_ham_im = join(path_dir_results, 'Ham_{}_im'.format(i))
        file_ham_re = join(path_dir_results, 'Ham_{}_re'.format(i))
        ham_im = 2.0 * css
        ham_re = np.diag(2.0 * energies)  # convert to Rydbergs
        write_pyxaid_format(ham_im, file_ham_im)
        write_pyxaid_format(ham_re, file_ham_re)
        ham_files.append((file_ham_im, file_ham_re))

    return ham_files


def change_mol_units(mol, factor=angs2au):
    """change the units of the molecular coordinates"""
    newMol = []
    for atom in mol:
        coord = list(map(lambda x: x * factor, atom.xyz))
        newMol.append(AtomXYZ(atom.symbol, coord))
    return newMol


def retrieve_hdf5_data(path_hdf5, paths_to_prop):
    """
    Read Numerical properties from ``paths_hdf5``.

    :params path_hdf5: Path to the hdf5 file
    :type path_hdf5: string
    :returns: numerical array
    """
    try:
        with h5py.File(path_hdf5, 'r') as f5:
            if isinstance(paths_to_prop, list):
                return  [f5[path][...] for path in paths_to_prop]
            else:
                return f5[paths_to_prop][...]
    except KeyError:
        msg = "There is not {} stored in the HDF5\n".format(paths_to_prop)
        with open('cp2k_out', 'a') as f:
            f.write(msg)
        raise KeyError(msg)
    except FileNotFoundError:
        msg = "there is not HDF5 file containing the numerical results"
        raise RuntimeError(msg)
    
