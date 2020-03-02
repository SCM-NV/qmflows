from typing import Type, Mapping

from scm import plams
plams.Job = plams.core.basejob.Job
plams.ORCAJob = plams.interfaces.thirdparty.orca.ORCAJob

from qmflows.packages.packages import Package

JOB_MAP: Mapping[Type[plams.Job], str] = {
    plams.Cp2kJob: 'generic2CP2K.json',
    plams.ADFJob: 'generic2ADF.json',
    plams.DFTBJob: 'generic2DFTB',
    plams.GamessJob: 'generic2gamess.json',
    plams.ORCAJob: 'generic2ORCA.sjon'
}


class PackageWrapper(Package):
    def __init__(self, job_type: Type[plams.Job], pkg_name: Optional[str] = None) -> None:
        """Initialize this instance."""
        if pkg_name is None:
            pkg_name = job_type.__class__.__name__.lower().rstrip('job')
        super(CP2K, self).__init__(pkg_name)
        self.generic_dict_file = JOB_MAP.get(job_type, 'generic2None.json')

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
                raise RuntimeError("cell parameter:{}\nformat not recognized")

            return s

        funs = {'cell_parameters': write_cell_parameters,
                'cell_angles': write_cell_angles}

        # Function that handles the special keyword
        f = funs.get(key)

        if f is not None:
            f(settings, value, mol, key)