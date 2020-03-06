"""CP2K input/output handling."""
__all__ = ['CP2K_Result', 'cp2k']

import os
from os.path import join
from typing import Optional, Union, Any, Dict, ClassVar
from warnings import warn

from scm import plams

from ..parsers.cp2KParser import parse_cp2k_warnings
from ..settings import Settings
from ..warnings_qmflows import cp2k_warnings
from .packages import (Package, Result, package_properties,
                       parse_output_warnings, WarnMap)

# ====================================<>=======================================
charge_dict: Dict[str, int] = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7,
    'Ne': 8, 'Na': 9, 'Mg': 10, 'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7,
    'Ar': 8, 'K': 9, 'Ca': 10, 'Sc': 11, 'Ti': 12, 'V': 13, 'Cr': 14, 'Mn': 15,
    'Fe': 16, 'Co': 17, 'Ni': 18, 'Cu': 11, 'Zn': 12, 'Ga': 13, 'Ge': 4, 'As': 5,
    'Se': 6, 'Br': 7, 'Kr': 8, 'Rb': 9, 'Sr': 10, 'Y': 11, 'Zr': 12, 'Nb': 13,
    'Mo': 14, 'Tc': 15, 'Ru': 16, 'Rh': 17, 'Pd': 18, 'Ag': 11, 'Cd': 12,
    'In': 13, 'Sn': 4, 'Sb': 5, 'Te': 6, 'I': 7, 'Xe': 8, 'Cs': 9, 'Ba': 10,
    'Hf': 12, 'Ta': 13, 'W': 14, 'Re': 15, 'Os': 16, 'Ir': 17, 'Pt': 18,
    'Au': 11, 'Hg': 12, 'Tl': 13, 'Pb': 4, 'Bi': 5, 'Po': 6, 'At': 7, 'Rn': 8}
# ======================================<>====================================
__all__ = ['cp2k']


class CP2K(Package):
    """This class setup the requirement to run a CP2K Job <https://www.cp2k.org/>.

    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    **cp2k** function should be called.

    """

    generic_dict_file: ClassVar[str] = 'generic2CP2K.yaml'

    def __init__(self) -> None:
        super().__init__("cp2k")

    @staticmethod
    def run_job(settings: Settings, mol: plams.Molecule,
                job_name: str = 'cp2k_job',
                work_dir: Union[None, str, os.PathLike] = None,
                **kwargs: Any) -> 'CP2K_Result':
        """Call the Cp2K binary using plams interface.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param mol: molecular Geometry
        :type mol: plams Molecule
        :param hdf5_file: Path to the HDF5 file that contains the numerical results.
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

        # Create a Plams job
        job = plams.interfaces.thirdparty.cp2k.Cp2kJob(
            name=job_name, settings=cp2k_settings, molecule=mol)
        r = job.run()

        work_dir = work_dir if work_dir is not None else job.path

        warnings = parse_output_warnings(job_name, r.job.path,
                                         parse_cp2k_warnings, cp2k_warnings)

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        result = CP2K_Result(cp2k_settings, mol, job_name, r.job.path, dill_path,
                             work_dir=work_dir, status=job.status, warnings=warnings)
        return result

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Create the settings input for complex cp2k keys.

        :param settings: Job Settings.
        :type settings: :class:`~qmflows.Settings`
        :param key: Special key declared in ``settings``.
        :param value: Value store in ``settings``.
        :param mol: molecular Geometry
        :type mol: plams Molecule

        """
        def write_cell_angles(s: Settings, value: Any,
                              mol: plams.Molecule, key: str) -> Settings:
            """The angles of the cell is a 3-dimensional list.

            &SUBSYS
              &CELL
                ABC [angstrom] 5.958 7.596 15.610
                ALPHA_BETA_GAMMA 81.250 86.560 89.800
              &END CELL

            """
            if value is not None:
                angles = '{} {} {}'.format(*value)
                s.specific.cp2k.force_eval.subsys.cell.alpha_beta_gamma = angles

            return s

        def write_cell_parameters(s: Settings, value: Any,
                                  mol: plams.Molecule, key: str) -> Settings:
            """The cell parameter can be a list of lists containing the ABC parameter.

            For example: ::

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
            def fun(xs):
                return '{} {} {}'.format(*xs)

            if not isinstance(value, list):
                abc = [value] * 3
                abc_cell = ' [angstrom] {} {} {}'.format(*abc)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc_cell
            elif isinstance(value[0], list):
                a, b, c = value
                s.specific.cp2k.force_eval.subsys.cell.A = fun(a)
                s.specific.cp2k.force_eval.subsys.cell.B = fun(b)
                s.specific.cp2k.force_eval.subsys.cell.C = fun(c)
            elif isinstance(value, list):
                abc = ' [angstrom] {} {} {}'.format(*value)
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
            else:
                raise RuntimeError(f"cell parameter:{value!r}\nformat not recognized")

            return s

        funs = {'cell_parameters': write_cell_parameters,
                'cell_angles': write_cell_angles}

        # Function that handles the special keyword
        f = funs.get(key)

        if f is not None:
            f(settings, value, mol, key)
        else:
            warn(f'Generic keyword {key!r} not implemented for package CP2K')


class CP2K_Result(Result):
    """Class providing access to CP2K result."""

    def __init__(self, settings: Optional[Settings],
                 molecule: Optional[plams.Molecule],
                 job_name: str,
                 dill_path: Union[None, str, os.PathLike] = None,
                 plams_dir: Union[None, str, os.PathLike] = None,
                 work_dir: Union[None, str, os.PathLike] = None,
                 status: str = 'successful',
                 warnings: Optional[WarnMap] = None) -> None:
        """Initialize this instance."""
        super().__init__(settings, molecule, job_name, plams_dir, dill_path,
                         work_dir=work_dir, properties=package_properties['cp2k'],
                         status=status, warnings=warnings)


# instance
cp2k = CP2K()
