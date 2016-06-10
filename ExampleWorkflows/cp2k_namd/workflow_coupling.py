__author__ = "Felipe Zapata"

# ================> Python Standard  and third-party <==========
from os.path import join
import os
import plams

# ==================> Internal modules <==========
from components import (calculate_mos, create_dict_CGFs, create_point_folder,
                         split_file_geometries)
from noodles import gather, schedule

from qmworks import run, Settings
from qmworks.parsers import parse_string_xyz
from nac.schedule.scheduleCoupling import (lazy_schedule_couplings,
                                           schedule_transf_matrix,
                                           write_hamiltonians)

# ==============================> Main <==================================


def generate_pyxaid_hamiltonians(package_name, project_name, all_geometries,
                                 cp2k_args, guess_args=None,
                                 calc_new_wf_guess_on_points=[0],
                                 path_hdf5=None, enumerate_from=0):
    """
    Use a md trajectory to generate the hamiltonian components to tun PYXAID
    nmad.

    :param package_name: Name of the package to run the QM simulations.
    :type  package_name: String
    :param project_name: Folder name where the computations
    are going to be stored.
    :type project_name: String
    :param all_geometries: List of string cotaining the molecular geometries
    numerical results.
    :type path_traj_xyz: [String]
    :param package_args: Specific settings for the package
    :type package_args: dict
    :param use_wf_guess_each: number of Computations that used a previous
    calculation as guess for the wave function.
    :type use_wf_guess_each: Int
    :param enumerate_from: Number from where to start enumerating the folders
    create for each point in the MD
    :type enumerate_from: Int
    :returns: None
    """
    #  Environmental Variables
    cwd = os.path.realpath(".")
    
    basisName = cp2k_args.basis
    work_dir = os.path.join(cwd, project_name)
    if path_hdf5 is None:
        path_hdf5 = os.path.join(work_dir, "quantum.hdf5")

    # Create Work_dir if it does not exist
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

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
    traj_folders = create_point_folder(work_dir, len(all_geometries),
                                       enumerate_from)

    # prepare Cp2k Jobs
    # Point calculations Using CP2K
    mo_paths_hdf5 = calculate_mos(package_name, all_geometries, work_dir,
                                  path_hdf5, traj_folders, cp2k_args,
                                  guess_args, calc_new_wf_guess_on_points,
                                  enumerate_from)

    # Calculate Non-Adiabatic Coupling
    # Number of Coupling points calculated with the MD trajectory
    nPoints = len(all_geometries) - 2
    promise_couplings = [calculate_coupling(i, path_hdf5, dictCGFs,
                                            all_geometries,
                                            mo_paths_hdf5, hdf5_trans_mtx,
                                            enumerate_from,
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
                                       path_dir_results=path_hamiltonians,
                                       enumerate_from=enumerate_from)

    hams_files = run(promise_files)

    print(hams_files)
# ==============================> Tasks <=====================================


def calculate_coupling(i, path_hdf5, dictCGFs, all_geometries, mo_paths,
                       hdf5_trans_mtx, enumerate_from, output_folder=None):
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
    :param enumerate_from: Number from where to start enumerating the folders
    create for each point in the MD
    :type enumerate_from: Int
    :returns: promise to path to the Coupling inside the HDF5
    """
    j, k = i + 1, i + 2
    geometries = all_geometries[i], all_geometries[j], all_geometries[k]

    return lazy_schedule_couplings(i, path_hdf5, dictCGFs, geometries, mo_paths,
                                   dt=1, hdf5_trans_mtx=hdf5_trans_mtx,
                                   output_folder=output_folder,
                                   enumerate_from=enumerate_from)
# ============<>===============


def main():
    """
    Initialize the arguments to compute the nonadiabatic coupling matrix for
    a given MD trajectory.
    """
    plams.init()
    project_name = 'ET_Pb79S44'

    # create Settings for the Cp2K Jobs
    cp2k_args = Settings()
    cp2k_args.basis = "DZVP-MOLOPT-SR-GTH"
    cp2k_args.potential = "GTH-PBE"
    cp2k_args.cell_parameters = [50.0] * 3
    cp2k_args.specific.cp2k.force_eval.dft.scf.added_mos = 100
    cp2k_args.specific.cp2k.force_eval.dft.scf.diagonalization.jacobi_threshold = 1e-6

    # Setting to calculate the WF use as guess
    cp2k_OT = Settings()
    cp2k_OT.basis = "DZVP-MOLOPT-SR-GTH"
    cp2k_OT.potential = "GTH-PBE"
    cp2k_OT.cell_parameters = [50.0] * 3
    cp2k_OT.specific.cp2k.force_eval.dft.scf.scf_guess = 'atomic'
    cp2k_OT.specific.cp2k.force_eval.dft.scf.ot.minimizer = 'DIIS'
    cp2k_OT.specific.cp2k.force_eval.dft.scf.ot.n_diis = 7
    cp2k_OT.specific.cp2k.force_eval.dft.scf.ot.preconditioner = 'FULL_SINGLE_INVERSE'
    cp2k_OT.specific.cp2k.force_eval.dft.scf.added_mos = 0
    cp2k_OT.specific.cp2k.force_eval.dft.scf.eps_scf = 5e-06

    # Path to the MD geometries
    path_traj_xyz = "./trajectory_4000-5000.xyz"

    # Work_dir
    scratch = "/scratch-shared"
    scratch_path = join(scratch, project_name)
    if not os.path.exists(scratch_path):
        os.makedirs(scratch_path)

    # HDF5 path
    path_hdf5 = join(scratch_path, 'quantum.hdf5')

    # all_geometries type :: [String]
    geometries = split_file_geometries(path_traj_xyz)

    # Named the points of the MD starting from this number
    enumerate_from = 0

    # Calculate new Guess in each Geometry
    pointsGuess = [enumerate_from + i for i in range(len(geometries))]

    # Hamiltonian computation
    generate_pyxaid_hamiltonians('cp2k', project_name, geometries, cp2k_args,
                                 guess_args=cp2k_OT,
                                 calc_new_wf_guess_on_points=pointsGuess,
                                 path_hdf5=path_hdf5,
                                 enumerate_from=enumerate_from)

    print("PATH TO HDF5:{}\n".format(path_hdf5))
    plams.finish()

if __name__ == "__main__":
    main()

