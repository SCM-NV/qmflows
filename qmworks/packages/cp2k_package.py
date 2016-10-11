

__all__ = ['cp2k']

# =======>  Standard and third party Python Libraries <======
from functools import partial
from warnings import warn
from os.path import join

import fnmatch
import h5py
import os
import plams

# ==================> Internal modules <====================
from qmworks.common import InputKey
from qmworks.hdf5 import cp2k2hdf5
from qmworks.packages.packages import Package, Result
from qmworks.parsers import read_cp2k_number_of_orbitals
from qmworks.settings import Settings
from qmworks.utils import lookup

# ====================================<>=======================================
charge_dict = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'S': 6, 'Cl': 7,
               'Se': 6, 'Cd': 12, 'Pb': 4, 'Br': 7, 'Cs': 9, 'Si': 4}
# ======================================<>====================================


class CP2K(Package):
    """
    This class setup the requirement to run a CP2K Job <https://www.cp2k.org/>.
    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    **cp2k** function should be called.
    """
    def __init__(self):
        super(CP2K, self).__init__("cp2k")
        self.generic_dict_file = 'generic2CP2K.json'

    def prerun(self):
        pass

    def run_job(self, settings, mol, work_dir=None, project_name=None,
                hdf5_file="quantum.hdf5", input_file_name=None,
                out_file_name=None, store_in_hdf5=True,
                nHOMOS=None, nLUMOS=None, job_name='cp2k_job'):
        """
        Call the Cp2K binary using plams interface.

        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the
        numerical results.
        :type hdf5_file: String
        :param input_file_name: Optional name for the input.
        :type input_file_name: String
        :param out_file_name: Optional name for the output.
        :type out_file_name: String
        :param store_in_hdf5: wether to store the output arrays in HDF5 format.
        :type store_in_hdf5: Bool
        """
        cp2k_settings = Settings()
        cp2k_settings.input = settings.specific.cp2k
        job = plams.Cp2kJob(name=job_name, settings=cp2k_settings,
                            molecule=mol)
        runner = plams.JobRunner(parallel=True)
        r = job.run(runner)
        r.wait()

        work_dir = work_dir if work_dir is not None else job.path
        output_file = join(job.path, job._filename('out'))

        if store_in_hdf5:
            dump_to_hdf5(hdf5_file, settings, work_dir, output_file, nHOMOS,
                         nLUMOS, project_name=project_name)

        return CP2K_Result(cp2k_settings, mol, job_name, r.job.path, work_dir,
                           path_hdf5=hdf5_file, project_name=project_name)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        """
        Create the settings input for complex cp2k keys

        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param key: Special key declared in ``settings``.
        :param value: Value store in ``settings``.
        :param mol: molecular Geometry
        :type mol: plams Molecule
        """

        def write_cell_parameters(s, value, mol, key):
            """
            The cell parameter can be a list of lists containing the
            ABC parameter like ::

            &SUBSYS
               &CELL
               A  16.11886919    0.07814137      -0.697284243
               B  -0.215317662   4.389405268     1.408951791
               C  -0.216126961   1.732808365     9.748961085
               PERIODIC XYZ
               &END
            .....

            The cell parameter can also be a scalar for ABC like ::

            &SUBSYS
            &CELL
            ABC [angstrom] 12.74 12.74 12.74
            PERIODIC NONE
            &END CELL
            """
            if not isinstance(value, list):
                abc = [value] * 3
                abc_cell = ' [angstrom] {} {} {}'.format(*abc)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc_cell
            elif isinstance(value[0], list):
                a, b, c = value  #
                fun = lambda xs: '{} {} {}'.format(*xs)
                s.specific.cp2k.force_eval.subsys.cell.A = fun(a)
                s.specific.cp2k.force_eval.subsys.cell.B = fun(b)
                s.specific.cp2k.force_eval.subsys.cell.C = fun(c)
                s.specific.cp2k.force_eval.subsys.cell.periodic = 'xyz'
            else:
                a, b, c = value  # Pattern match list
                abc = ' [angstrom] {} {} {}'.format(a, b, c)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
                s.specific.cp2k.force_eval.subsys.cell.periodic = 'xyz'

            return s

        def expand_basis_set(s, prefix, mol, key):
            """
            CP2k has a sspecial format for the basis set, For more
            information have a look at
            `basis <https://www.cp2k.org/basis_sets?s[]=basis>`.
            For a Molecule that contains only carbon and oxygen atoms,
            the basis set declaration is given by,
            >>> &FORCE_EVAL
                    .......
                    &SUBSYS
                        &KIND  C
                            BASIS_SET  DZVP-MOLOPT-SR-GTH-q4
                            POTENTIAL  GTH-PBE-q4
                        &END C
                        &KIND  H
                            BASIS_SET  DZVP-MOLOPT-SR-GTH-q1
                            POTENTIAL  GTH-PBE-q1
                        &END H
                    &END SUBSYS
                & END FORCE_EVALXS
            Where DZVP-MOLOPT-SR-GTH is the name of the basis and q4, q1
            correspond to the charge associated with that atom
            (e.g. 4 for carbon, 1 for hydrogen).
            """

            def symbols2charge(s):
                q = charge_dict[s]
                return 'q{}'.format(q)

            symbols = set([at.symbol for at in mol.atoms])
            qs = list(map(symbols2charge, symbols))
            for symb, q in zip(symbols, qs):
                name = '{}-{}'.format(prefix, q)
                if key == 'basis':
                    s.specific.cp2k.force_eval.subsys.kind[symb]["basis_set"] = name
                elif key == 'potential':
                    s.specific.cp2k.force_eval.subsys.kind[symb]["potential"] = name
            return s

        funs = {'basis': expand_basis_set, 'potential': expand_basis_set,
                'cell_parameters': write_cell_parameters}

        # Function that handles the special keyword
        f = lookup(funs, key)

        if f is not None:
            return f(settings, value, mol, key)
        else:
            warn(UserWarning('Keyword ' + key + ' doesn\'t exist'))


class CP2K_Result(Result):
    """
    Class providing access to CP2K result.
    """
    def __init__(self, settings, molecule, job_name, plams_dir, work_dir=None,
                 path_hdf5=None, project_name=None,
                 properties='data/dictionaries/propertiesCP2K.json'):
        super().__init__(settings, molecule, job_name, plams_dir,
                         work_dir=work_dir, path_hdf5=path_hdf5,
                         project_name=project_name, properties=properties)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        """
        Create a :class:`~CP2K_Result` instance using the data serialized in
        a dictionary.

        :param cls:
        :param settings: Job Settings.
        :param molecule: molecular Geometry.
        :param job_name: Name of the job.
        :param plams_dir: Absolute path to plams output folder
        :param archive: dictionary containing the paths to the input/output
        folders.
        :param path_hdf5: Path to the HDF5 file that contains the numerical
        results.
        """
        fun = partial(lookup, archive)
        plams_dir, work_dir, path_hdf5 = list(map(fun, ["plams_dir", "work_dir",
                                                        "path_hdf5"]))
        return CP2K_Result(settings, molecule, job_name, plams_dir, work_dir,
                           path_hdf5, project_name)

    def get_property(self, prop, section=None):
        pass

    def __getattr__(self, prop):
        """Returns a section of the results.
        The property is extracted from a  result file, which is recursively
        search for in the CP2K settings

        Example:
        ..
            overlap_matrix = result.overlap
        """
        relative_cwd = self.archive['work_dir'].split('/')[-1]
        hdf5_path_to_prop = join(self.project_name, relative_cwd)
        sections = self.prop_dict[prop]
        paths_to_prop = list(map(lambda x: join(hdf5_path_to_prop, x),
                                 sections))

        return paths_to_prop


cp2k = CP2K()


def dump_to_hdf5(file_h5, settings, work_dir, output_file, nHOMOS,
                 nLUMOS, project_name=None):
    """
    Store the result in HDF5 format.

    :param file_h5: Path to the HDF5 file that contains the
    numerical results.
    :type file_h5: String
    :param settings: Job Settings.
    :type settings: :class:`~qmworks.Settings`
    :param work_dir: Path to the folders where the calculation is carried out.
    :tpye work_dir: String
    :param output_file: Absolute path to plams output file.
    """
    def match_file(pattern):
        """ Cp2k append the suffix .Log to the output files """
        s = "*{}*Log".format(pattern)
        xs = list(filter(lambda x: fnmatch.fnmatch(x, s),
                         os.listdir(work_dir)))
        if xs:
            return xs[0]
        else:
            return None

    def get_value_recursively(st, xs):
        """
        :param xs: List of keys
        :type xs: String List
        """
        s = st.copy()
        for x in xs:
            s = s.get(x)
            if s is None:
                break
        return s

    def get_file_path(xs):
        """
        Search for a result file requested in the settings.
        CP2K renames thew files appending a number an a `Log` to the
        end of the filename.
        """
        path = get_value_recursively(settings, xs)
        if path is None or os.path.exists(path):
            return path
        else:  # The software renamed the filename given by the user
            root, file_pattern = os.path.split(path)
            real_name = match_file(file_pattern)
            return join(root, real_name)

    settings_file_MO = ["specific", "cp2k", "force_eval", "dft",
                        "print", "mo", "filename"]
    settings_file_overlap = ["specific", "cp2k", "force_eval", "dft",
                             "print", "ao_matrices", "filename"]

    # Arguments to store the properties in HDF5
    nOccupied, nOrbitals, nOrbFuns = read_cp2k_number_of_orbitals(output_file)
    path_MO = get_file_path(settings_file_MO)
    path_overlap = get_file_path(settings_file_overlap)

    # Paths inside the HDF5 file
    keys = []
    files_to_remove = []
    if path_MO is not None:
        relative_cwd = work_dir.split('/')[-1]
        pathEs = join(project_name, relative_cwd, "cp2k/mo/eigenvalues")
        pathCs = join(project_name, relative_cwd, "cp2k/mo/coefficients")
        k = InputKey('orbitals',
                     [path_MO, nOrbitals, nOrbFuns, pathEs, pathCs,
                      nOccupied, nHOMOS, nLUMOS])
        keys.append(k)
        # Remove this file after it has been processed
        files_to_remove.append(path_MO)

    if path_overlap is not None:
        path_mtx_overlap = join(work_dir, "cp2k/overlap")
        keys.append(InputKey('overlap',
                             [path_overlap, nOrbitals, path_mtx_overlap]))
        files_to_remove.append(path_overlap)

    # Calling the qmworks-HDF5 API
    with h5py.File(file_h5, chunks=True) as f5:
        cp2k2hdf5(f5, keys)

    # Remove the text output files
    for x in files_to_remove:
        os.remove(x)
