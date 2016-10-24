
# ================> Python Standard  and third-party <==========
from collections import namedtuple
from nose.plugins.attrib import attr
from noodles import schedule  # Workflow Engine
from os.path import join
from qmworks import (run, Settings, templates)
from qmworks.packages import cp2k
from qmworks.utils import (chunksOf)

import fnmatch
import h5py
import os
import plams
import shutil
# ===================================<>========================================
JobFiles = namedtuple("JobFiles", ("get_xyz", "get_inp", "get_out", "get_MO"))


@attr('slow')
def test_ethylene():
    """
    run a single point calculation using CP2K and store the MOs.
    """
    home = os.path.expanduser('~')  # HOME Path
    scratch_path = join(home, '.test_qmworks')
    if not os.path.exists(scratch_path):
        os.makedirs(scratch_path)
    try:
        fun_ethylene(scratch_path)
    finally:
        # remove tmp data and clean global config
        shutil.rmtree(scratch_path)


def fun_ethylene(scratch_path):
    """
    Test Ethylene singlw
    """
    project_name = 'ethylene'

    # create Settings for the Cp2K Jobs
    s = Settings()
    s.basis = "DZVP-MOLOPT-SR-GTH"
    s.potential = "GTH-PBE"
    s.cell_parameters = [12.74] * 3
    dft = s.specific.cp2k.force_eval.dft
    dft.scf.added_mos = 20
    dft.scf.eps_scf = 1e-4

    dft['print']['ao_matrices']['overlap'] = ''
    dft['print']['ao_matrices']['filename'] = join(scratch_path, 'overlap.out')

    # Copy the basis and potential to a tmp file
    shutil.copy('test/test_files/BASIS_MOLOPT', scratch_path)
    shutil.copy('test/test_files/GTH_POTENTIALS', scratch_path)
    # Cp2k configuration files
    basiscp2k = join(scratch_path, 'BASIS_MOLOPT')
    potcp2k = join(scratch_path, 'GTH_POTENTIALS')
    cp2k_config = {"basis": basiscp2k, "potential": potcp2k}

    # HDF5 path
    path_hdf5 = join(scratch_path, 'ethylene.hdf5')

    # all_geometries type :: [String]
    path_xyz = 'test/test_files/ethylene.xyz'
    geometries = split_file_geometries(path_xyz)

    # Input/Output Files
    file_xyz = join(scratch_path, 'coordinates.xyz')
    file_inp = join(scratch_path, 'ethylene.inp')
    file_out = join(scratch_path, 'ethylene.out')
    file_MO = join(scratch_path, 'mo_coeffs.out')

    files = JobFiles(file_xyz, file_inp, file_out, file_MO)

    schedule_job = schedule(prepare_job_cp2k)

    promise = schedule_job(geometries[0], files, s, scratch_path,
                           project_name=project_name, hdf5_file=path_hdf5,
                           wfn_restart_job=None, store_in_hdf5=True, nHOMOS=25,
                           nLUMOS=25, package_config=cp2k_config)

    cp2k_result = run(promise)

    path_properties = cp2k_result.orbitals

    with h5py.File(path_hdf5) as f5:
        assert(all(p in f5 for p in path_properties))


def prepare_job_cp2k(geometry, files, settings, work_dir,
                     project_name=None, hdf5_file=None, wfn_restart_job=None,
                     store_in_hdf5=True, nHOMOS=20, nLUMOS=20,
                     package_config=None):
    """
    Fills in the parameters for running a single job in CP2K.

    :param geometry: Molecular geometry stored as String
    :type geometry: String
    :param files: Tuple containing the IO files to run the calculations
    :type files: nameTuple
    :parameter settings: Dictionary contaning the data to
    fill in the template
    :type      settings: Dict
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
    job_settings = prepare_cp2k_settings(geometry, files, settings,
                                         work_dir, wfn_restart_job,
                                         store_in_hdf5, package_config)
    print("CP2K Settings: ", job_settings)

    return cp2k(job_settings, plams.Molecule(files.get_xyz), work_dir=work_dir,
                project_name=project_name, hdf5_file=hdf5_file,
                input_file_name=files.get_inp,
                out_file_name=files.get_out, store_in_hdf5=store_in_hdf5,
                nHOMOS=nHOMOS, nLUMOS=nLUMOS)


def prepare_cp2k_settings(geometry, files, cp2k_args, work_dir,
                          wfn_restart_job, store_in_hdf5, cp2k_config):
    """
    Fills in the parameters for running a single job in CP2K.

    :param geometry: Molecular geometry stored as String
    :type geometry: String
    :param files: Tuple containing the IO files to run the calculations
    :type files: nameTuple
    :parameter settings: Dictionary contaning the data to
    fill in the template
    :type  settings: ~qmworks.Settings
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
    cp2k_args.specific.cp2k['global']['project'] = 'ethylene'

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
    return list(map(''.join, chunksOf(xss, numat + 2)))
