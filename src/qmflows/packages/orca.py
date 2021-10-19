"""Orca input/output bookkeeping."""
__all__ = ['orca']

import os
from os.path import join
from warnings import warn
from typing import Any, Union, Optional, ClassVar, List, Type

import numpy as np
from scm import plams

from .packages import Package, Result, load_properties
from ..parsers.orca_parser import parse_molecule
from ..settings import Settings
from ..type_hints import Final, _Settings
from ..utils import get_tmpfile_name
from ..warnings_qmflows import Key_Warning

# ============================= Orca ==========================================


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results."""

    prop_mapping: ClassVar[_Settings] = load_properties('ORCA', prefix='properties')

    @property
    def molecule(self) -> Optional[plams.Molecule]:
        """Retrieve the molecule from the output file."""
        if self.status in {'crashed', 'failed'}:
            return None

        plams_dir = self.archive["plams_dir"]
        if plams_dir is None:
            return None
        else:
            file_name = join(plams_dir, f'{self.job_name}.out')
            return parse_molecule(file_name, self._molecule)


class ORCA(Package):
    """A class for preparing the input for running a Orca job using both PLAMS and templates.

    It also does the manangement of the input/output files resulting
    from running Orca and returns a Results object that containing the methods
    and data required to retrieve the output.

    """

    generic_mapping: ClassVar[_Settings] = load_properties('ORCA', prefix='generic2')
    result_type: ClassVar[Type[ORCA_Result]] = ORCA_Result

    def __init__(self, pkg_name: str = "orca") -> None:
        super().__init__(pkg_name)

    @classmethod
    def run_job(cls, settings: Settings, mol: plams.Molecule,
                job_name: str = "ORCAjob",
                work_dir: "None | str | os.PathLike[str]" = None,
                validate_output: bool = True,
                **kwargs: Any) -> ORCA_Result:

        orca_settings = Settings()
        orca_settings.input = settings.specific.orca

        # Running Orca with Plams
        job = plams.interfaces.thirdparty.orca.ORCAJob(
            molecule=mol, settings=orca_settings, name=job_name)
        result = job.run()

        # Relative job path
        relative_plams_path = join(*str(result.job.path).split(os.sep)[-2:])

        # Absolute path to the .dill file
        dill_path = join(job.path, f'{job.name}.dill')

        return cls.result_type(orca_settings, mol, result.job.name, dill_path,
                               plams_dir=relative_plams_path, status=job.status)

    @staticmethod
    def handle_special_keywords(settings: Settings, key: str,
                                value: Any, mol: plams.Molecule) -> None:
        """Translate generic keywords to their corresponding Orca keywords."""
        def inithess(value: Any) -> Settings:
            """Generate an seperate file containing the initial Hessian matrix.

            It is used as guess for the computation.

            """
            # Convert Hessian to numpy array
            value = value if isinstance(value, np.ndarray) else np.array(value)

            def format_atom(atom: plams.Atom) -> str:
                symbol, mass, coords = atom.symbol, atom.mass, atom.coords
                return '{:2s}{:12.4f}{:14.6f}{:14.6f}{:14.6f}\n'.format(symbol, mass, *coords)

            def format_hessian(dim, hess: Union[List[List[float]], np.ndarray]) -> str:
                """Format numpy array to Orca matrix format."""
                ret = ''
                for i in range((dim - 1) // 6 + 1):
                    n_columns = min(6, dim - 6 * i)
                    ret += '         '
                    ret += ' '.join('{:10d}'.format(v + 6 * i)
                                    for v in range(n_columns))
                    ret += '\n'
                    for j in range(dim):
                        ret += '{:7d}     '.format(j)
                        ret += ' '.join('{:10.6f}'.format(hess[v + 6 * i][j])
                                        for v in range(n_columns))
                        ret += '\n'
                return ret

            # Check Hessian dimension
            if value.ndim == 1:
                dim = int(value.size ** 0.5)
                hess = np.reshape(value, (dim, dim))
            else:
                dim = value.shape[0]
                hess = value

            # Header
            hess_str = '\n$orca_hessian_file\n\n$hessian\n' + str(dim) + '\n'
            # Actual Hessian
            hess_str += format_hessian(dim, hess)
            # Atoms header
            hess_str += '\n$atoms\n'
            # Atoms coordinates
            hess_str += str(len(mol)) + '\n'
            hess_str += ''.join(format_atom(atom) for atom in mol.atoms)
            # The end
            hess_str += '\n\n$end\n'

            # Store the hessian in the plams_dir
            hess_path = get_tmpfile_name("ORCA_hessian_")
            with open(hess_path, "w") as hess_file:
                hess_file.write(hess_str)

            settings.specific.orca.geom.InHess = "read"
            settings.specific.orca.geom.InHessName = f'{hess_path.as_posix()}'

            return settings

        def constraint(value: Any) -> None:
            cons = ''
            if isinstance(value, Settings):
                for k, v in value.items():
                    ks = k.split()
                    atoms = [int(a) - 1 for a in ks[1:]]
                    if ks[0] == 'dist' and len(ks) == 3:
                        cons += '{{ B {:d} {:d} {:.2f} C }}'.format(*atoms, v)
                    elif ks[0] == 'angle' and len(ks) == 4:
                        cons += '{{ A {:d} {:d} {:d} {:.2f} C }}'.format(
                            *atoms, v)
                    elif ks[0] == 'dihed' and len(ks) == 5:
                        cons += '{{ D {:d} {:d} {:d} {:d} {:.2f} C }}'.format(
                            *atoms, v)
                    else:
                        warn(
                            f'Invalid constraint key: {k}', category=Key_Warning)
            settings.specific.orca.geom.Constraints._end = cons

        def freeze(value: List[Union[int, str]]) -> None:
            if not isinstance(value, list):
                msg = f'selected_atoms {value} is not a list'
                raise RuntimeError(msg)
            cons = ''
            if isinstance(value[0], int):
                for a in value:
                    cons += '{{ C {:d} C }}'.format(a - 1)
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol in value:
                        cons += '{{ C {:d} C }}'.format(a)
            settings.specific.orca.geom.Constraints._end = cons

        def selected_atoms(value: List[Union[int, str]]) -> None:
            if not isinstance(value, list):
                raise RuntimeError(f'selected_atoms {value} is not a list')
            cons = ''
            if isinstance(value[0], int):
                for a in range(len(mol)):
                    if a + 1 not in value:
                        cons += '{{ C {:d} C }}'.format(a)
            else:
                for a in range(len(mol)):
                    if mol[a + 1].symbol not in value:
                        cons += '{{ C {:d} C }}'.format(a)
            settings.specific.orca.geom.Constraints._end = cons

        # Available translations
        functions = {'inithess': inithess,
                     'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'constraint': constraint}
        if key in functions:
            functions[key](value)
        else:
            warn(f'Generic keyword {key!r} not implemented for package ORCA',
                 category=Key_Warning)


#: An instance :class:`ORCA`.
#: Only one instance of this class should exist at any given momemt;
#: *i.e.* this value is a singleton.
orca: Final[ORCA] = ORCA()
