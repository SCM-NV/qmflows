"""CP2K input/output handling."""
__all__ = ['CP2K_Result', 'cp2k']

import os
from os.path import join
from typing import Any, Union, Iterable, ClassVar, Type
from warnings import warn

import numpy as np
from scm import plams

from .packages import Package, Result, parse_output_warnings, load_properties
from ..parsers.cp2KParser import parse_cp2k_warnings
from ..settings import Settings
from ..warnings_qmflows import cp2k_warnings, Key_Warning
from ..type_hints import Final, _Settings

__all__ = ['cp2k']


class CP2K_Result(Result):
    """Class providing access to CP2K result."""

    prop_mapping: ClassVar[_Settings] = load_properties('CP2K', prefix='properties')

    @property
    def molecule(self) -> "None | plams.Molecule":
        """Return the current geometry.

        If the job is an optimization, try to read the ` *-pos-1.xyz` file.
        Otherwise return the input molecule.
        """
        try:
            return self.get_property('geometry')
        except FileNotFoundError:
            return self._molecule


class CP2K(Package):
    """A Package subclass for running `CP2K Jobs <https://www.cp2k.org/>`_.

    It uses plams together with the templates to generate the stucture input
    and also uses Plams to invoke the binary CP2K code.
    This class is not intended to be called directly by the user, instead the
    :data:`cp2k` function should be called.

    """

    generic_mapping: ClassVar[_Settings] = load_properties('CP2K', prefix='generic2')
    result_type: ClassVar[Type[CP2K_Result]] = CP2K_Result

    def __init__(self, pkg_name: str = "cp2k") -> None:
        super().__init__(pkg_name)

    @classmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule,
                job_name: str = 'cp2k_job',
                work_dir: "None | str | os.PathLike[str]" = None,
                validate_output: bool = True,
                **kwargs: Any) -> CP2K_Result:
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
        cp2k_settings.executable = settings.get("executable", "cp2k.popt")

        # Create a Plams job
        job = plams.interfaces.thirdparty.cp2k.Cp2kJob(
            name=job_name, settings=cp2k_settings, molecule=mol)
        r = job.run()

        work_dir = work_dir if work_dir is not None else job.path

        if validate_output:
            warnings = parse_output_warnings(
                job_name, r.job.path, parse_cp2k_warnings, cp2k_warnings
            )
        else:
            warnings = None

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        result = cls.result_type(cp2k_settings, mol, job_name, dill_path=dill_path,
                                 plams_dir=r.job.path, work_dir=work_dir,
                                 status=job.status, warnings=warnings)
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
        def write_cell_angles(s: Settings, value: list,
                              mol: plams.Molecule, key: str) -> None:
            """Write the Angles for the cell.

            The angles of the cell is a 3-dimensional list.
            &SUBSYS
              &CELL
                ABC [angstrom] 5.958 7.596 15.610
                ALPHA_BETA_GAMMA 81.250 86.560 89.800
              &END CELL

            """
            if value is not None:
                angles = '{} {} {}'.format(*value)
                s.specific.cp2k.force_eval.subsys.cell.ALPHA_BETA_GAMMA = angles

        def write_cell_parameters(s: Settings, value: Any,
                                  mol: plams.Molecule, key: str) -> None:
            """Write the cell parameters.

            The cell parameter can be a list of lists containing the ABC parameter.

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
            def fun(xs: Iterable[Any]) -> str:
                return '{:} {:} {:}'.format(*xs)

            ar = np.asarray(value, dtype=np.float64)
            if ar.ndim == 0:
                abc = f' [angstrom] {fun(np.repeat(ar, 3))}'
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
            elif ar.ndim == 1:
                abc = f' [angstrom] {fun(ar)}'
                s.specific.cp2k.force_eval.subsys.cell.ABC = abc
            elif ar.ndim == 2:
                a, b, c = ar
                s.specific.cp2k.force_eval.subsys.cell.A = fun(a)
                s.specific.cp2k.force_eval.subsys.cell.B = fun(b)
                s.specific.cp2k.force_eval.subsys.cell.C = fun(c)
            else:
                raise RuntimeError(
                    f"cell parameter:{value!r}\nformat not recognized")

        def write_periodic(s: Settings, value: Any, mol: plams.Molecule, key: str) -> None:
            """Set the keyword for periodic calculations."""
            s.specific.cp2k.force_eval.subsys.cell.periodic = value

        def ignore_keyword(s: Settings, value: Any, mol: plams.Molecule, key: str) -> None:
            pass

        funs = {'cell_parameters': write_cell_parameters,
                'cell_angles': write_cell_angles,
                'periodic': write_periodic,
                'executable': ignore_keyword}

        # Function that handles the special keyword
        f = funs.get(key)

        if f is not None:
            f(settings, value, mol, key)
        else:
            warn(f'Generic keyword {key!r} not implemented for package CP2K',
                 category=Key_Warning)


#: An instance :class:`CP2K`.
#: Only one instance of this class should exist at any given momemt;
#: *i.e.* this value is a singleton.
cp2k: Final[CP2K] = CP2K()
