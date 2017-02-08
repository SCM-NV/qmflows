# =======>  Standard and third party Python Libraries <======
from warnings import warn
import plams

# ==================> Internal modules <====================
from qmworks.packages.packages import (Package, package_properties, Result)
from qmworks.settings import Settings

# ====================================<>=======================================
charge_dict = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8, 'Na': 9, 'Mg': 10, 'Al': 3,
               'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8, 'K': 9, 'Ca': 10, 'Sc' : 11, 'Ti': 12, 'V': 13, 'Cr': 14, 'Mn': 15,
               'Fe': 16, 'Co': 17, 'Ni': 18, 'Cu': 11, 'Zn': 12, 'Ga': 13, 'Ge': 4, 'As': 5, 'Se': 6, 'Br': 7, 'Kr': 8, 'Rb': 9,
               'Sr': 10, 'Y': 11, 'Zr': 12, 'Nb': 13, 'Mo': 14, 'Tc': 15, 'Ru': 16, 'Rh': 17, 'Pd': 18, 'Ag': 11, 'Cd': 12,
               'In': 13, 'Sn': 4, 'Sb': 5, 'Te': 6, 'I': 7, 'Xe': 8, 'Cs': 9, 'Ba': 10, 'Hf': 12, 'Ta': 13, 'W': 14, 'Re': 15,
               'Os': 16, 'Ir': 17, 'Pt': 18, 'Au': 11, 'Hg': 12, 'Tl': 13, 'Pb': 4, 'Bi': 5, 'Po': 6, 'At': 7, 'Rn': 8}
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

        # Add molecular coordinates
        m = format_coord_xyz(mol) + '{:>8}'.format('&END')
        cp2k_settings.input.force_eval.subsys['&COORD'] = m

        # Create a Plams job
        job = plams.interfaces.cp2k.Cp2kJob(name=job_name,
                                            settings=cp2k_settings,
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

        def write_cell_angles(s, value, mol, key):
            """
            The angles of the cell is a 3-dimensional list ::

            &SUBSYS
              &CELL
                ABC [angstrom] 5.958 7.596 15.610
                ALPHA_BETA_GAMMA 81.250 86.560 89.800
              &END CELL
            """
            angles = '{} {} {}'.format(*value)
            s.specific.cp2k.force_eval.subsys.cell.alpha_beta_gamma = angles

            return s

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
            elif isinstance(value, list):
                abc = ' [angstrom] {} {} {}'.format(*value)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
            elif isinstance(value[0], list):
                a, b, c = value  #
                fun = lambda xs: '{} {} {}'.format(*xs)
                s.specific.cp2k.force_eval.subsys.cell.A = fun(a)
                s.specific.cp2k.force_eval.subsys.cell.B = fun(b)
                s.specific.cp2k.force_eval.subsys.cell.C = fun(c)
            else:
                msg = "cell parameter:{}\nformat not recognized"
                RuntimeError(msg)

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
                'cell_parameters': write_cell_parameters,
                'cell_angles': write_cell_angles}

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
        return CP2K_Result(settings, molecule, job_name, plams_dir.path,
                           work_dir=work_dir,
                           properties=package_properties['cp2k'],
                           status=status)


def format_coord_xyz(mol):
    """
    """
    xs  = ''.join('{}  {: 12.8e}  {: 12.8e}  {: 12.8e}\n'.format(at.symbol, *at.coords)
                  for at in mol)

    return '\n' + xs

cp2k = CP2K()
