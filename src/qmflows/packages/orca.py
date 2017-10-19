__all__ = ['orca']

from os.path import join
from qmflows.settings import Settings
from qmflows.packages.packages import (Package, package_properties, Result, get_tmpfile_name)
from qmflows.parsers.orca_parser import parse_molecule
from scm import plams
from warnings import warn

import builtins
import numpy as np
# ============================= Orca ==========================================


class ORCA(Package):
    """
    This class prepare the input to run a Orca job using both Plams and
    templates. It also does the manangement of the input/output files resulting
    from running Orca and returns a Results object that containing the methods
    and data required to retrieve the output.
    """
    def __init__(self):
        super().__init__("orca")
        self.generic_dict_file = 'generic2ORCA.json'

    def prerun(self):
        pass

    @staticmethod
    def run_job(settings, mol, job_name="ORCAjob"):

        orca_settings = Settings()
        orca_settings.input = settings.specific.orca

        # Running Orca with Plams
        job = plams.interfaces.orca.ORCAJob(molecule=mol,
                                            settings=orca_settings,
                                            name=job_name)
        result = job.run()

        # Relative job path
        relative_plams_path = '/'.join(result.job.path.split('/')[-2:])

        return ORCA_Result(orca_settings, mol, result.job.name,
                           plams_dir=relative_plams_path, status=job.status)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Translate generic keywords to their corresponding Orca keywords.
        """

        def inithess(value):
            """
            Generate an seperate file containing the initial Hessian matrix used as
            guess for the computation.
            """
            # Convert Hessian to numpy array
            value = value if isinstance(value, np.ndarray) else np.array(value)

            def format_atom(atom):
                symbol, mass, coords = atom.symbol, atom._getmass(), atom.coords
                return '{:2s}{:12.4f}{:14.6f}{:14.6f}{:14.6f}\n'.format(symbol, mass, *coords)

            def format_hessian(dim, hess):
                """ Format numpy array to Orca matrix format """
                ret = ''
                for i in range((dim - 1) // 6 + 1):
                    n_columns = min(6, dim - 6 * i)
                    ret += '         '
                    ret += ' '.join('{:10d}'.format(v + 6 * i) for v in range(n_columns))
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
            hess_path = get_tmpfile_name()
            with open(hess_path, "w") as hess_file:
                hess_file.write(hess_str)

            settings.specific.orca.geom.InHess = "read"
            settings.specific.orca.geom.InHessName = '"' + hess_path + '"'

            return settings

        def constraint(value):
            cons = ''
            if isinstance(value, Settings):
                for k, v in value.items():
                    ks = k.split()
                    atoms = [int(a) - 1 for a in ks[1:]]
                    if ks[0] == 'dist' and len(ks) == 3:
                        cons += '{{ B {:d} {:d} {:f} C }}'.format(*atoms, v)
                    elif ks[0] == 'angle' and len(ks) == 4:
                        cons += '{{ A {:d} {:d} {:d} {:f} C }}'.format(*atoms, v)
                    elif ks[0] == 'dihed' and len(ks) == 5:
                        cons += '{{ D {:d} {:d} {:d} {:d} {:f} C }}'.format(*atoms, v)
                    else:
                        warn('Invalid constraint key: ' + k)
            settings.specific.orca.geom.Constraints._end = cons

        def freeze(value):
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            cons = ''
            if isinstance(value[0], int):
                for a in value:
                    cons += '{{ C {:d} C }}'.format(a - 1)
            else:
                for a in range(len(mol)):
                    if mol[a+1].symbol in value:
                        cons += '{{ C {:d} C }}'.format(a)
            settings.specific.orca.geom.Constraints._end = cons

        def selected_atoms(value):
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            cons = ''
            if isinstance(value[0], int):
                for a in range(len(mol)):
                    if a + 1 not in value:
                        cons += '{{ C {:d} C }}'.format(a)
            else:
                for a in range(len(mol)):
                    if mol[a+1].symbol not in value:
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
            msg = 'Keyword ' + key + ' not implemented for package ORCA'
            warn(msg)


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results"""

    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 status='done', warnings=None):
        properties = package_properties['orca']
        super().__init__(settings, molecule, job_name=job_name,
                         plams_dir=plams_dir, properties=properties,
                         status=status, warnings=warnings)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status, warnings):
        plams_dir = archive["plams_dir"].path
        return ORCA_Result(settings, molecule, job_name, plams_dir, status, warnings)

    @property
    def molecule(self):
        """ Retrieve the molecule from the output file"""
        if self.status not in ['crashed', 'failed']:
            plams_dir = self.archive["plams_dir"].path
            file_name = join(plams_dir, '{}.out'.format(self.job_name))
            return parse_molecule(file_name, self._molecule)
        else:
            return None

orca = ORCA()
