__author__ = "Felipe Zapata"

# ================> Python Standard  and third-party <==========
from collections import namedtuple
from itertools import starmap
from os.path import join

import h5py
import itertools
import os
import plams

# ==================> Internal modules <==========
from nac.basisSet.basisNormalization import createNormalizedCGFs
from noodles import gather, schedule

from qmworks.fileFunctions import search_environ_var
from qmworks import run
from qmworks.common import InputKey
from qmworks.hdf5.quantumHDF5 import cp2k2hdf5
from qmworks.packages import cp2k_farming
from qmworks.parsers import parse_string_xyz
from qmworks.utils import chunksOf, flatten
from scheduleCp2k import (prepare_cp2k_settings, prepare_farming_cp2k_settings)
from scheduleCoupling import (lazy_schedule_couplings, schedule_transf_matrix,
                              write_hamiltonians)

# ==============================<>=========================
# Tuple contanining file paths
JobFiles = namedtuple("JobFiles", ("get_xyz", "get_inp", "get_out", "get_MO"))

# ==============================> Main <=======================================


def main():
    plams.init()
    project_name = 'NAC'
    package_args = {'added_mos': 100, 'basisName': "DZVP-MOLOPT-SR-GTH",
                    'potname': "GTH-PBE",
                    'abc_cell': [[16.11886919, 0.07814137, -0.697284243],
                                 [-0.215317662, 4.389405268, 1.408951791],
                                 [-0.216126961, 1.732808365, 9.748961085]]}

    # Path to the MD geometries
    geometries = "./data/traj_6_points.xyz"

    # Hamiltonian computation
    calculate_pyxaid_hams('cp2k', project_name, geometries, package_args)
    plams.finish()


def calculate_pyxaid_hams(package_name, project_name, path_traj_xyz,
                          package_args):
    """
    Use a md trajectory to generate the hamiltonian components to tun PYXAID
    nmad.

    :param package_name: Name of the package to run the QM simulations.
    :type  package_name: String
    :param project_name: Folder name where the computations
    are going to be stored.
    :type project_name: String
    :param path_traj_xyz: Path to the md trajectory stored as XYZ.
    numerical results.
    :type path_traj_xyz: String
    :param package_args: Specific settings for the package
    :type package_args: dict
    :returns: None
    """
    #  Environmental Variables
    cwd = os.path.realpath(".")
    
    basisName = package_args['basisName']
    work_dir = os.path.join(cwd, project_name)
    path_hdf5 = os.path.join(work_dir, "quantum.hdf5")

    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # all_geometries type :: [String]
    all_geometries = split_file_geometries(path_traj_xyz)

    # Generate a list of tuples containing the atomic label
    # and the coordinates to generate
    # the primitive CGFs
    atoms = parse_string_xyz(all_geometries[0])
    dictCGFs = create_dict_CGFs(path_hdf5, basisName, atoms)

    # Calculcate the matrix to transform from cartesian to spherical
    # representation of the overlap matrix
    hdf5_trans_mtx = schedule_transf_matrix(path_hdf5, atoms,
                                            basisName, work_dir,
                                            packageName=package_name)
    
    # Create a folder for each point the the dynamics
    traj_folders = creat_point_folder(work_dir, len(all_geometries))

    # prepare Cp2k Jobs
    # Point calculations Using CP2K
    mo_paths_hdf5 = calculate_mos(package_name, all_geometries, work_dir,
                                  path_hdf5, traj_folders, package_args)

    # Calculate Non-Adiabatic Coupling
    # Number of Coupling points calculated with the MD trajectory
    nPoints = len(all_geometries) - 2
    promise_couplings = [calculate_coupling(i, path_hdf5, dictCGFs,
                                            all_geometries,
                                            mo_paths_hdf5, hdf5_trans_mtx,
                                            output_folder=work_dir)
                         for i in range(nPoints)]
    path_couplings = gather(*promise_couplings)

    # Write the results in PYXAID format
    path_hamiltonians = join(work_dir, 'hamiltonians')
    if not os.path.exists(path_hamiltonians):
        os.makedirs(path_hamiltonians)

    # Inplace scheduling of write_hamiltonians function.
    # Equivalent to add @schedule on top of the function
    schedule_write_ham = schedule(write_hamiltonians)
        
    promise_files = schedule_write_ham(path_hdf5, work_dir, mo_paths_hdf5,
                                       path_couplings, nPoints,
                                       path_dir_results=path_hamiltonians)
    hams_files = run(promise_files)

    print(hams_files)

# ==============================> Tasks <=====================================


def calculate_coupling(i, path_hdf5, dictCGFs, all_geometries,
                       mo_paths, hdf5_trans_mtx, output_folder=None):
    """
    Calculate the non-adiabatic coupling using 3 consecutive set of MOs in
    a dynamics. Explicitly declares that each Coupling Depends in
    three set of MOs.

    :param i: nth coupling calculation.
    :type i: Int
    :param path_hdf5: Path to the HDF5 file that contains the
    numerical results.
    :type path_hdf5: String
    :paramter dictCGFS: Dictionary from Atomic Label to basis set
    :type     dictCGFS: Dict String [CGF],
              CGF = ([Primitives], AngularMomentum),
              Primitive = (Coefficient, Exponent)
    :param all_geometries: list of molecular geometries
    :type all_geometries: String list
    :param mo_paths: Path to the MO coefficients and energies in the
    HDF5 file.
    :type mo_paths: [String]
    :param hdf5_trans_mtx: path to the transformation matrix in the HDF5 file.
    :type hdf5_trans_mtx: String
    """
    j, k = i + 1, i + 2
    geometries = all_geometries[i], all_geometries[j], all_geometries[k]

    return lazy_schedule_couplings(i, path_hdf5, dictCGFs, geometries, mo_paths,
                                   dt=1, hdf5_trans_mtx=hdf5_trans_mtx,
                                   output_folder=output_folder)


def calculate_mos(package_name, all_geometries, work_dir, path_hdf5, folders,
                  package_args):
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
    :returns: path to nodes in the HDF5 file to MO energies and MO coefficients.
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

        if not os.path.exists(path_hdf5):
            return None
        else:
            with h5py.File(path_hdf5, 'r') as f5:
                if isinstance(paths_to_prop, list):
                    pred = all(path in f5 for path in paths_to_prop)
                else:
                    pred = paths_to_prop in f5

                return paths_to_prop if pred else  None

    def prepare_jobs(job_files, start_index=0):
        """
        Prepare the CP2K job setting for the individual jobs that
        are going to be run in Farming mode. The 2nd job uses the
        WF of the first job as a guess, the 3rd the WF from 2nd
        and so on.
        """
        gss = all_geometries[start_index:]
        files = jobs_files[start_index:]
        work_dirs = folders[start_index:]
        indexes = itertools.count(start_index)
        farming = False
        restart = None

        # The first job does not have a restart to read from
        jobs = []
        for gs, fs, wd, k in zip(gss, files, work_dirs, indexes):
            s = prepare_cp2k_settings(gs, fs, package_args, k, wd,
                                      hdf5_file=path_hdf5,
                                      wfn_restart_job=restart,
                                      farming_use_guess=farming)
            restart = join(wd, 'point_{}-RESTART.wfn'.format(k))
            farming = True
            jobs.append(s)

        return jobs
    
    # Create all the job file names associated with a point in the MD
    jobs_files = list(starmap(create_file_names, enumerate(folders)))
    input_files = [f.get_inp for f in jobs_files]
    output_files = [f.get_out for f in jobs_files]
    xyz_files = [f.get_xyz for f in jobs_files]

    # Check if the MOs orbitals are already available at the HDF5 file
    paths_to_prop = [search_data_in_hdf5(k) for k in range(len(folders))]

    if all(x is not None for x in paths_to_prop):
        return paths_to_prop
    else:
        # Create the setting for the pending jobs
        jobs = prepare_jobs(jobs_files)
    
        farming_settings = prepare_farming_cp2k_settings(folders, input_files)
    
        r = cp2k_farming(farming_settings, None, hdf5_file=path_hdf5,
                         work_dirs=folders, input_files=input_files,
                         output_files=output_files, coordinates_files=xyz_files,
                         jobs_to_harvest=jobs)
        return r.orbitals


def creat_point_folder(work_dir, n):
    """
    Create a new folder for each point in the MD trajectory.

    :returns: Paths lists.
    """
    folders = []
    for k in range(n):
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


def create_file_names(i, work_dir):
    """
    Creates a namedTuple with the name of the three files used
    for each point in the trajectory
    
    :returns: Namedtuple containing the IO files
    """
    file_xyz = join(work_dir, 'coordinates_{}.xyz'.format(i))
    file_inp = join(work_dir, 'point_{}.inp'.format(i))
    file_out = join(work_dir, 'point_{}.out'.format(i))
    file_MO = join(work_dir, 'mo_coeff_{}.out'.format(i))

    return JobFiles(file_xyz, file_inp, file_out, file_MO)


# FIXME: Extended to other QM packages
def create_dict_CGFs(path_hdf5, basisname, xyz):
    """
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


def splitAtNone(xs):
    """Yield Elements of the list until one of the elements is None"""
    for i, x in enumerate(xs):
        if x is None:
            break

    return xs[:i]


# ============<>===============
if __name__ == "__main__":
    main()
