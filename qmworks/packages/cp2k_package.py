

__all__ = ['cp2k', 'cp2k_farming', 'CP2K', 'CP2K_Farming', 'CP2K_Result']

# =======>  Standard and third party Python Libraries <======
from noodles import files
from os.path import join

import fnmatch
import h5py
import os
import pkg_resources as pkg
import plams
# ==================> Internal modules <====================
from qmworks.common import InputKey
from qmworks.fileFunctions import json2Settings
from qmworks.hdf5 import cp2k2hdf5
from qmworks.packages.packages import Package, Result
from qmworks.parsers import read_cp2k_number_of_orbitals
from qmworks.settings import Settings

# ====================================<>========================================
charge_dict = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'S': 6, 'Cl': 7,
               'Se': 6, 'Cd': 12, 'Pb': 4}

# ======================================<>======================================


class CP2K(Package):
    """
    This class setup the requirement to run a CP2K Job <https://www.cp2k.org/>.
    It uses plams together with the templates to generate the stucture input and
    also uses Plams to invoke the binary CP2K code.
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
                nHOMOS=100, nLUMOS=100):
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
        job = plams.Cp2kJob(settings=cp2k_settings, molecule=mol)
        runner = plams.JobRunner(parallel=True)
        r = job.run(runner)
        r.wait()

        work_dir = work_dir if work_dir is not None else job.path
        output_file = join(job.path, job._filename('out'))

        if store_in_hdf5:
            self.dump_to_hdf5(hdf5_file, settings, work_dir, output_file, nHOMOS,
                              nLUMOS)

        return CP2K_Result(cp2k_settings, mol, r.job.path, work_dir, hdf5_file,
                           project_name)

    def postrun(self):
        pass

    def dump_to_hdf5(self, file_h5, settings, work_dir, output_file, nHOMOS,
                     nLUMOS):
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
        files = os.listdir(work_dir)

        def match_file(pattern):
            """ Cp2k append the suffix .Log to the output files """
            s = "*{}*Log".format(pattern)
            xs = list(filter(lambda x: fnmatch.fnmatch(x, s), files))
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
                else:
                    s
            return s

        def get_file_path(xs):
            """
            Search for a result file requested in the settings.
            CP2K renames thew files appending a number an a `Log` to the
            end of the filename.
            """

            path = get_value_recursively(settings, xs)
            if path is None:
                return None
            else:
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
            pathEs = join(work_dir, "cp2k/mo/eigenvalues")
            pathCs = join(work_dir, "cp2k/mo/coefficients")
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
                s.specific.cp2k.force_eval.subsys.cell.PERIODIC = 'XYZ'
            else:
                a, b, c = value  # Pattern match list
                abc = ' [angstrom] {} {} {}'.format(a, b, c)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
                s.specific.cp2k.force_eval.subsys.cell.PERIODIC = 'XYZ'

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

        return funs[key](settings, value, mol, key)


class CP2K_Farming(CP2K):
    """
    Run CP2K Job in Groups according to:
    <https://manual.cp2k.org/trunk/CP2K_INPUT/FARMING.html>
    """
    def run_job(self, settings, mol=None, hdf5_file="quantum.hdf5",
                work_dirs=None, input_files=None, output_files=None,
                coordinates_files=None, jobs_to_harvest=None, nGroups=None,
                initial_guess=None):
        """
        Runs a CP2K Farming Job

        :param settings: Farming Job Settings
        :type settings: :class:`~qmworks.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the
        numerical results.
        :type hdf5_file: String
        :param work_dirs: List of directories containing the information to
        run a single job (e.g. coordinates).
        :type work_dirs: String List
        :param input_file_names: Names of the input files to be included in the
        farming job.
        :type input_file_names: String List
        :param out_file_names: Output File Names for the different jobs
        included in the farming.
        :type out_file_names: String List
        :param nGroups: Number of Jobs to run in parallel inside the Farming.
        :type nGroups: Int
        :returns: `~CP2K_Result`
        """
        # Execute initial_guess
        if initial_guess is not None:
                    initial_guess.orbitals
        
        # Create the input files for the Jobs to be run in Farming mode
        mols = [plams.Molecule(xyz) for xyz in coordinates_files]
        job_settings = [self.generic2specific(j, m)
                        for j, m in zip(jobs_to_harvest, mols)]

        # Use plams to generate the input of the jobs to farm
        plams_settings = [Settings() for x in job_settings]
        for cp2k_plams, s, input_file in zip(plams_settings, job_settings,
                                             input_files):
            cp2k_plams.input = s.specific.cp2k
            job = plams.Cp2kJob(settings=cp2k_plams)

            # Generate Input with Plams to use it during the farming
            input = job.get_input()
            with open(input_file, 'w') as f:
                f.write(input)

        # pass arguments in Plams input settings format
        cp2k_settings = Settings()
        cp2k_settings.input = settings.specific.cp2k
        farm_job = plams.Cp2kJob(name='farming', settings=cp2k_settings)

        # Create a Farming Job
        runner = plams.JobRunner(parallel=True)
        r = farm_job.run(runner)
        r.wait()

        # Name of the folder were the output is stored
        folder_names = [f.split('/')[-1] for f in work_dirs]

        # Dump the numerical results to HDF5
        for s, outFile, folder in zip(job_settings, output_files, work_dirs):
            self.dump_to_hdf5(hdf5_file, s, folder, outFile)

        return CP2K_Farming_Result(cp2k_settings, mol, r.job.path, work_dirs,
                                   hdf5_file)

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
        def expand_farming_jobs(s, values, mol, key):
            """
            The Input for a Farming CP2K job resemble the following structure,

             &GLOBAL
                PROJECT farming
                PROGRAM FARMING
                RUN_TYPE NONE
             &END GLOBAL

             &FARMING
               GROUP_SIZE 1

               &JOB
                 DIRECTORY dir-1
                 INPUT_FILE_NAME job1.inp
                 JOB_ID 1
               &END JOB

               &JOB
                 DEPENDENCIES 1
                 DIRECTORY dir-2
                 INPUT_FILE_NAME job2.inp
                 JOB_ID 2
               &END JOB

               ...........................
               ...........................

               &JOB
                 DEPENDENCIES 31
                 DIRECTORY dir-32
                 INPUT_FILE_NAME job32.inp
                 JOB_ID 32
               &END JOB

             &END FARMING
            """
            job_ids, job_names, workDirs = [list(t) for t in zip(*values)]
            s.specific.cp2k.farming.job.directories = workDirs
            s.specific.cp2k.farming.job.input_file_names = job_names
            s.specific.cp2k.farming.job.job_ids = job_ids

            return s
                    
        funs = {'farming_jobs': expand_farming_jobs}

        try:
            return super(CP2K_Farming,
                         self).handle_special_keywords(settings, key, value, mol)
        except KeyError:
            return funs[key](settings, value, mol, key)


class CP2K_Result(Result):
    """
    Class providing access to CP2K result.
  
    :param settings:
    """
    def __init__(self, settings, molecule, plams_dir, work_dir, file_h5,
                 project_name):
        """
        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        """
        self.settings = settings
        self._molecule = molecule
        self.hdf5_file = file_h5
        properties = 'data/dictionaries/propertiesCP2K.json'
        xs = pkg.resource_string("qmworks", properties)
        self.prop_dict = json2Settings(xs)
        # self.archive = result_path
        self.archive = {"plams_dir": files.Path(plams_dir),
                        'work_dir': work_dir}
        self.project_name = project_name
        
    def as_dict(self):
        return {
            "settings": self.settings,
            "molecule": self._molecule,
            "filename": self.archive}

    @classmethod
    def from_dict(cls, settings, molecule, plams_dir=None, work_dir=None,
                  file_h5='quantum.hdf5'):
        """
        Create a :class:`~CP2K_Result` instance using the data serialized in
        a dictionary.

        :param cls:
        :param settings: Job Settings
        :param molecule: molecular Geometry
        :param plams_dir: Absolute path to plams output folder
        :param work_dir: Absolute path to the folder where the calculation
        was performed.
        :param file_h5: Path to the HDF5 file that contains the numerical results
        """
        return CP2K_Result(settings, molecule, plams_dir, work_dir, file_h5)

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
        sections = self.prop_dict[prop]
        # paths_to_prop = list(map(lambda x: join(self.archive['work_dir'], x),
        #                          sections))
        paths_to_prop = list(map(lambda x: join(self.project_name, x), sections))

        return paths_to_prop


class CP2K_Farming_Result(CP2K_Result):
    """
    """
    def __init__(self, settings, mol, plams_dir, work_dirs, file_h5):
        """
        :param settings: Job Settings.
        :type settings: :class:`~qmworks.Settings`
        :param mol: molecular Geometry.
        :type mol: plams Molecule
        :param plams_dir: Absolute path to plams output folder
        :type plams_dir: String
        :param work_dir: Absolute path to the folder where the calculation
        was performed.
        :type work_dirs: String
        :param file_h5: Path to the HDF5 file that contains the numerical results.
        :type file_h5: String
        """
        super(CP2K_Farming_Result, self).__init__(settings, mol, plams_dir,
                                                  work_dirs, file_h5)
        self.archive = {"plams_dir": files.Path(plams_dir),
                        'work_dirs': work_dirs}
        
    def __getattr__(self, prop):
        """Returns a section of the results.
        The property is extracted from a  result file, which is recursively
        search for in the CP2K settings

        Example:
        ..
            overlap_matrix = result.overlap
        """
        sections = self.prop_dict[prop]
        work_dirs = self.archive['work_dirs']
        paths_to_prop = list(map(lambda folder:
                                 list(map(lambda x: join(folder, x), sections)), work_dirs))

        return paths_to_prop


cp2k = CP2K()
cp2k_farming = CP2K_Farming()
