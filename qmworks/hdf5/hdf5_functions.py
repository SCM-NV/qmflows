from .quantumHDF5 import StoreasHDF5
from os.path import join

# ===============>


def dump_to_hdf5(data, package_name, file_h5, project_name='workflow',
                 job_name=None, property_to_dump='orbitals'):
    """
    Store the result in HDF5 format.

    :param file_h5: Path to the HDF5 file that contains the
    numerical results.
    :type file_h5: String
    :param work_dir: Path to the folders where the calculation is carried out.
    :tpye work_dir: String
    :param output_file: Absolute path to plams output file.
    :returns: None
    """
    job_name = job_name if job_name is not None else "job"
    store_hdf5 = StoreasHDF5(file_h5, package_name)
    if property_to_dump == 'orbitals':
        es = join(package_name, "mo", "eigenvalues")
        css = join(package_name, "mo", "coefficients")
        pathEs = join(project_name, job_name, es)
        pathCs = join(project_name, job_name, css)

    for p, d in zip([pathEs, pathCs], [data.eigenVals, data.coeffs]):
        store_hdf5.funHDF5(p, d)
