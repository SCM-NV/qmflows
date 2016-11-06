# =======>  Standard and third party Python Libraries <======
from warnings import warn
import plams

# ==================> Internal modules <====================
from qmworks.packages.packages import (Package, package_properties, Result)
from qmworks.settings import Settings

# ====================================<>=======================================
charge_dict = {'H': 1, 'C': 4, 'N': 5, 'O': 6, 'S': 6, 'Cl': 7,
               'Se': 6, 'Cd': 12, 'Pb': 4, 'Br': 7, 'Cs': 9, 'Si': 4}
# ======================================<>====================================
__all__ = ['cp2k']


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

    @staticmethod
    def run_job(settings, mol, job_name='cp2k_job',
                work_dir=None, **kwargs):
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
        # Yet another work directory

        # Input modifications
        cp2k_settings = Settings()
        cp2k_settings.input = settings.specific.cp2k
        job = plams.Cp2kJob(name=job_name, settings=cp2k_settings,
                            molecule=mol)
        r = job.run()

        work_dir = work_dir if work_dir is not None else job.path

        result = CP2K_Result(cp2k_settings, mol, job_name, r.job.path,
                             work_dir, status=job.status)

        return result

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
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
        f = funs.get(key)

        if f is not None:
            return f(settings, value, mol, key)
        else:
            msg = 'Keyword ' + key + ' doesn\'t exist'
            warn(msg)


class CP2K_Result(Result):
    """
    Class providing access to CP2K result.
    """
    def __init__(self, settings, molecule, job_name, plams_dir, work_dir=None,
                 properties=package_properties['cp2k'],
                 status='successful'):
        super().__init__(settings, molecule, job_name, plams_dir,
                         work_dir=work_dir, properties=properties,
                         status=status)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status):
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
        plams_dir, work_dir = list(map(archive.get, ["plams_dir", "work_dir"]))
        return CP2K_Result(settings, molecule, job_name, plams_dir, work_dir,
                           status)

cp2k = CP2K()
