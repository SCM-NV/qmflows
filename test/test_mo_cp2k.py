
# ================> Python Standard  and third-party <==========
from collections import namedtuple
from noodles import schedule  # Workflow Engine
from os.path import join
from qmworks import (Settings, templates)
from qmworks.packages import cp2k

import fnmatch
import getpass
import os
import plams

# ===================================<>========================================
JobFiles = namedtuple("JobFiles", ("get_xyz", "get_inp", "get_out", "get_MO"))


def test_ethylene():
    """
    run a single point calculation using CP2K and store the MOs.
    """
    plams.init()
    project_name = 'ethylene'

    # create Settings for the Cp2K Jobs
    s = Settings()
    s.basis = "DZVP-MOLOPT-SR-GTH"
    s.potential = "GTH-PBE"
    s.cell_parameters = [12.74] * 3
    s.specific.cp2k.force_eval.dft.scf.added_mos = 50
    s.specific.cp2k.force_eval.dft.scf.diagonalization.jacobi_threshold = 1e-6

    # User variables
    home = os.path.expanduser('~')  # HOME Path
    username = getpass.getuser()
    # Work_dir
    scratch = "./"
    scratch_path = join(scratch, username, project_name)
    if not os.path.exists(scratch_path):
        os.makedirs(scratch_path)

    # Cp2k configuration files
    basiscp2k = join(home, "test/test_files/BASIS_MOLOPT")
    potcp2k = join(home, "test/test_files/GTH_POTENTIALS")
    cp2k_config = {"basis": basiscp2k, "potential": potcp2k}

    # HDF5 path
    path_hdf5 = join(scratch_path, 'quantum.hdf5')

    # all_geometries type :: [String]
    geometries = split_file_geometries(path_traj_xyz)


    # Input/Output Files
    file_xyz = join(scratch_path, 'coordinates.xyz')
    file_inp = join(scratch_path, 'ethylene.inp')
    file_out = join(scratch_path, 'ethylene.out')
    file_MO =  join(scratch_path, 'mo_coeffs.out')

    files = JobFiles(file_xyz, file_inp, file_out, file_MO)


    plams.finish()


def prepare_cp2k_settings(geometry, files, cp2k_args, k, work_dir,
                          wfn_restart_job, store_in_hdf5, cp2k_config):
    """
    Fills in the parameters for running a single job in CP2K.

    :param geometry: Molecular geometry stored as String
    :type geometry: String
    :param files: Tuple containing the IO files to run the calculations
    :type files: nameTuple
    :parameter dict_input: Dictionary contaning the data to
    fill in the template
    :type  dict_input: Dict
    :param k: nth Job
    :type k: Int
    :parameter work_dir: Name of the Working folder
    :type      work_dir: String
    :param wfn_restart_job: Path to *.wfn cp2k file use as restart file.
    :type wfn_restart_job: String
    :param store_in_hdf5: Wether or not the numerical result are store in HDF5.
    numerical results.
    :type store_in_hdf5: Bool
    :param cp2k_config:  Parameters required by cp2k.
    :type cp2k_config: Dict
   :returns: ~qmworks.Settings
    """
    # Search for the environmental variable BASISCP2K containing the path
    # to the Basis set folder

    basis_file = cp2k_config["basis"]
    potential_file = cp2k_config["potential"]

    force = cp2k_args.specific.cp2k.force_eval
    force.dft.basis_set_file_name = basis_file
    force.dft.potential_file_name = potential_file
    force.dft['print']['mo']['filename'] = files.get_MO
    force.subsys.topology.coord_file_name = files.get_xyz
    cp2k_args.specific.cp2k['global']['project'] = 'point_{}'.format(k)

    if wfn_restart_job is not None:
        output_dir = getattr(wfn_restart_job.archive['plams_dir'], 'path')
        xs = os.listdir(output_dir)
        wfn_file = list(filter(lambda x: fnmatch.fnmatch(x, '*wfn'), xs))[0]
        file_path = join(output_dir, wfn_file)
        force.dft.wfn_restart_file_name = file_path

    with open(files.get_xyz, 'w') as f:
        f.write(geometry)

    input_args = templates.singlepoint.overlay(cp2k_args)

    # Do not print the MO if they are not going to be stored in HDF5
    if not store_in_hdf5:
        del(input_args.specific.cp2k.force_eval.dft['print'])

    return input_args


@schedule
def prepare_job_cp2k(geometry, files, dict_input, k, work_dir,
                     project_name=None, hdf5_file=None, wfn_restart_job=None,
                     store_in_hdf5=True, nHOMOS=25, nLUMOS=25,
                     package_config=None):
    """
    Fills in the parameters for running a single job in CP2K.

    :param geometry: Molecular geometry stored as String
    :type geometry: String
    :param files: Tuple containing the IO files to run the calculations
    :type files: nameTuple
    :parameter dict_input: Dictionary contaning the data to
    fill in the template
    :type      dict_input: Dict
    :param k: nth Job
    :type k: Int
    :parameter work_dir: Name of the Working folder
    :type      work_dir: String
    :param hdf5_file: Path to the HDF5 file that contains the
    numerical results.
    :type hdf5_file: String
    :param wfn_restart_job: Path to *.wfn cp2k file use as restart file.
    :type wfn_restart_job: String
    :param farming_use_guess: Use a guess for the WF using a previous job.
    :type farming_use_guess: Bool
    :param nHOMOS: number of HOMOS to store in HDF5.
    :type nHOMOS: Int
    :param nLUMOS: number of HOMOS to store in HDF5.
    :type nLUMOS: Int
    :returns: ~qmworks.CP2K
    """
    job_settings = prepare_cp2k_settings(geometry, files, dict_input, k,
                                         work_dir, wfn_restart_job,
                                         store_in_hdf5, package_config)
    project_name = project_name if project_name is not None else work_dir

    return cp2k(job_settings, plams.Molecule(files.get_xyz), work_dir=work_dir,
                project_name=project_name, hdf5_file=hdf5_file,
                input_file_name=files.get_inp,
                out_file_name=files.get_out, store_in_hdf5=store_in_hdf5,
                nHOMOS=nHOMOS, nLUMOS=nLUMOS)
