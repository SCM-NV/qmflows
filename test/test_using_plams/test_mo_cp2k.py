
# ================> Python Standard  and third-party <==========
from collections import namedtuple
from nose.plugins.attrib import attr
from os.path import join
from plams import Molecule
from qmworks import (run, Settings, templates)
from qmworks.hdf5 import dump_to_hdf5
from qmworks.packages import cp2k

import h5py
import os
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
    Test Ethylene single
    """
    geometry = Molecule('test/test_files/ethylene.xyz')
    job_settings = prepare_cp2k_settings(geometry, scratch_path)

    cp2k_result = run(cp2k(job_settings, geometry, work_dir=scratch_path))

    # Path to the HDF5 file
    hdf5_file = 'quantum.hdf5'

    # Path to the molecular orbitals energies and  coefficients
    path_es = 'ethylene/cp2k_job/cp2k/mo/eigenvalues'
    path_css = 'ethylene/cp2k_job/cp2k/mo/coefficients'

    with h5py.File(hdf5_file) as f5:
        dump_to_hdf5(cp2k_result.orbitals, 'cp2k', f5,
                     project_name='ethylene',
                     job_name=cp2k_result.job_name,
                     property_to_dump='orbitals')
        energies = f5[path_es].value
        coefficients = f5[path_css].value

    print("Energy array shape: ", energies.shape)
    print("Coefficients array shape: ", coefficients.shape)
    assert (energies.shape == (40,)) and (coefficients.shape == (46, 40))


def prepare_cp2k_settings(geometry, work_dir):
    """
    Fills in the parameters for running a single job in CP2K.

    :param geometry: Molecular geometry stored as String
    :type geometry: plams.Molecule
    :param files: Tuple containing the IO files to run the calculations
    :type files: nameTuple
    :parameter settings: Dictionary contaning the data to
    fill in the template
    :type  settings: ~qmworks.Settings
    :parameter work_dir: Name of the Working folder
    :type      work_dir: String
    :param wfn_restart_job: Path to *.wfn cp2k file use as restart file.
    :type wfn_restart_job: String
    :param cp2k_config:  Parameters required by cp2k.
    :type cp2k_config: Dict
   :returns: ~qmworks.Settings
    """
    # Input/Output Files
    file_xyz = join(work_dir, 'coordinates.xyz')
    file_MO = join(work_dir, 'mo_coeffs.out')

    # create Settings for the Cp2K Jobs
    cp2k_args = Settings()
    cp2k_args.basis = "DZVP-MOLOPT-SR-GTH"
    cp2k_args.potential = "GTH-PBE"
    cp2k_args.cell_parameters = [12.74] * 3
    dft = cp2k_args.specific.cp2k.force_eval.dft
    dft.scf.added_mos = 20
    dft.scf.eps_scf = 1e-3
    dft["print"]["mo"]["mo_index_range"] = "7 46"
    dft.scf.diagonalization.jacobi_threshold = 1e-5

    # Copy the basis and potential to a tmp file
    shutil.copy('test/test_files/BASIS_MOLOPT', work_dir)
    shutil.copy('test/test_files/GTH_POTENTIALS', work_dir)
    # Cp2k configuration files

    force = cp2k_args.specific.cp2k.force_eval
    force.dft.basis_set_file_name = join(work_dir, 'BASIS_MOLOPT')
    force.dft.potential_file_name = join(work_dir, 'GTH_POTENTIALS')
    force.dft['print']['mo']['filename'] = file_MO
    force.subsys.topology.coord_file_name = file_xyz
    cp2k_args.specific.cp2k['global']['project'] = 'ethylene'

    # Write Molecular gemoetry
    geometry.write(file_xyz)

    return templates.singlepoint.overlay(cp2k_args)
