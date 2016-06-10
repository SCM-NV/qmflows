__author__ = "Felipe Zapata"

# ================> Python Standard  and third-party <==========
from collections import namedtuple
from noodles import gather
from os.path import join

import h5py
import os

# ==================> Internal modules <==========
from nac.basisSet.basisNormalization import createNormalizedCGFs
from nac.schedule.scheduleCp2k import prepare_job_cp2k
from qmworks.fileFunctions import search_environ_var
from qmworks.common import InputKey
from qmworks.hdf5.quantumHDF5 import cp2k2hdf5
from qmworks.utils import chunksOf, flatten


# ==============================<>=========================
# Tuple contanining file paths
JobFiles = namedtuple("JobFiles", ("get_xyz", "get_inp", "get_out", "get_MO"))

# ==============================> Tasks <=====================================


def calculate_mos(package_name, all_geometries, work_dir, path_hdf5, folders,
                  package_args, guess_args=None,
                  calc_new_wf_guess_on_points=[0], enumerate_from=0,
                  nHOMOS=100, nLUMOS=100):
    """
    Look for the MO in the HDF5 file if they do not exists calculate them by
    splitting the jobs in batches given by the ``restart_chunk`` variables.
    Only the first job is calculated from scratch while the rest of the
    batch uses as guess the wave function of the first calculation in
    the batch.

    :param all_geometries: list of molecular geometries
    :type all_geometries: String list
    :param work_dir: Path to the work directory
    :type work_dir: String
    :param path_hdf5: Path to the HDF5 file that contains the
    numerical results.
    :type path_hdf5: String
    :param folders: path to the directories containing the MO outputs
    :type folders: String list
    :param calc_new_wf_guess_on_points: Calculate a new Wave function guess in
    each of the geometries indicated. By Default only an initial guess is
    computed.
    :type calc_new_wf_guess_on_points: [Int]
    :param enumerate_from: Number from where to start enumerating the folders
    create for each point in the MD
    :type enumerate_from: Int
    :returns: path to nodes in the HDF5 file to MO energies
    and MO coefficients.
    """
    def create_properties_path(i):
        """
        Path inside HDF5 where the data is stored
        """
        rs = join(work_dir, 'point_{}'.format(i), package_name, 'mo')
        return [join(rs, 'eigenvalues'), join(rs, 'coefficients')]

    def search_data_in_hdf5(i):
        """
        Search if the node exists in the HDF5 file.
        """
        paths_to_prop = create_properties_path(i)

        with h5py.File(path_hdf5, 'r') as f5:
            if isinstance(paths_to_prop, list):
                pred = all(path in f5 for path in paths_to_prop)
            else:
                pred = paths_to_prop in f5

        return paths_to_prop if pred else None

    path_to_orbitals = []  # list to the nodes in the HDF5 containing the MOs

    # First calculation has no initial guess
    guess_job = None

    # calculate the rest of the job using the previous point as initial guess
    for j, gs in enumerate(all_geometries):
        k = j + enumerate_from
        point_dir = folders[j]
        job_files = create_file_names(point_dir, k)
        if k in calc_new_wf_guess_on_points:
            # Calculating initial guess
            guess_job = call_schedule_qm(package_name, guess_args, path_hdf5,
                                         point_dir, job_files, k, gs, nHOMOS,
                                         nLUMOS, guess_job=guess_job,
                                         store_in_hdf5=False)

        paths_to_prop = search_data_in_hdf5(k)

        # If the MOs are already store in the HDF5 format return the path
        # to them and skip the calculation
        if paths_to_prop is not None:
            path_to_orbitals.append(paths_to_prop)
        else:
            promise_qm = call_schedule_qm(package_name, package_args,
                                          path_hdf5, point_dir, job_files,
                                          k, gs, nHOMOS, nLUMOS,
                                          guess_job=guess_job)
            path_to_orbitals.append(promise_qm.orbitals)
            guess_job = promise_qm

    return gather(*path_to_orbitals)


def call_schedule_qm(packageName, package_args, path_hdf5, point_dir,
                     job_files, k, geometry, nHOMOS, nLUMOS, guess_job=None,
                     store_in_hdf5=True):
    """
    Call an external computational chemistry software to do some calculations

    :param package_name: Name of the package to run the QM simulations.
    :type  package_name: String
    :param package_args: Specific settings for the package
    :type package_args: Settings
    :param path_hdf5: Path to the HDF5 file that contains the
    numerical results.
    type path_hdf5: String
    :param point_dir: path to the directory where the output is written.
    :type point_dir: String
    :param job_files: Tuple containing the absolute path to IO files.
    :type job_files: NamedTuple
    :param k: current point being calculate in the MD
    :type k: Int
    :param geometry: Molecular geometry
    :type geometry: String
    :param nHOMOS: number of HOMOS to store in HDF5.
    :type nHOMOS: Int
    :param nLUMOS: number of HOMOS to store in HDF5.
    :type nLUMOS: Int
    :returns: promise QMWORK
    """
    prepare_and_schedule = {'cp2k': prepare_job_cp2k}

    job = prepare_and_schedule[packageName](geometry, job_files, package_args,
                                            k, point_dir, hdf5_file=path_hdf5,
                                            wfn_restart_job=guess_job,
                                            store_in_hdf5=store_in_hdf5,
                                            nHOMOS=nHOMOS, nLUMOS=nLUMOS)

    return job


def create_point_folder(work_dir, n, enumerate_from):
    """
    Create a new folder for each point in the MD trajectory.

    :returns: Paths lists.
    """
    folders = []
    for k in range(enumerate_from, n + enumerate_from):
        new_dir = join(work_dir, 'point_{}'.format(k))
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        folders.append(new_dir)

    return folders


def split_file_geometries(pathXYZ):
    """
    Reads a set of molecular geometries in xyz format and returns
    a list of string, where is element a molecular geometry
    
    :returns: String list containing the molecular geometries.
    """
    # Read Cartesian Coordinates
    with open(pathXYZ) as f:
        xss = f.readlines()

    numat = int(xss[0].split()[0])
    return list(map(flatten, chunksOf(xss, numat + 2)))


def create_dict_CGFs(path_hdf5, basisname, xyz):
    """
    FIXME: Extended to other QM packages
    If the Cp2k Basis are already stored in the hdf5 file continue,
    otherwise read and store them in the hdf5 file.
    """
    # Try to read the basis otherwise read it from a file
    with h5py.File(path_hdf5, chunks=True) as f5:
        try:
            f5["cp2k/basis"]
        except KeyError:
            # Search Path to the file containing the basis set
            pathBasis = search_environ_var('BASISCP2K')
            keyBasis = InputKey("basis", [pathBasis])
            cp2k2hdf5(f5, [keyBasis])             # Store the basis sets
        # Read the basis Set from HDF5 and calculate the CGF for each atom
        dictCGFs = createNormalizedCGFs(f5, basisname, 'cp2k', xyz)

    return dictCGFs


def create_file_names(work_dir, i):
    """
    Creates a namedTuple with the name of the 4 files used
    for each point in the trajectory
    
    :returns: Namedtuple containing the IO files
    """
    file_xyz = join(work_dir, 'coordinates_{}.xyz'.format(i))
    file_inp = join(work_dir, 'point_{}.inp'.format(i))
    file_out = join(work_dir, 'point_{}.out'.format(i))
    file_MO = join(work_dir, 'mo_coeff_{}.out'.format(i))

    return JobFiles(file_xyz, file_inp, file_out, file_MO)

# ============<>===============


