# =======>  Standard and third party Python Libraries <======
from os.path import join
from qmworks.settings import Settings
from qmworks.packages.packages import (Package, package_properties, Result)
from qmworks.parsers.orca_parser import parse_molecule
from warnings import warn

import builtins
import plams
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

        job = plams.ORCAJob(molecule=mol, settings=orca_settings,
                            name=job_name)
        result = job.run()

        return ORCA_Result(orca_settings, mol, result.job.name,
                           plams_dir=result.job.path, status=job.status)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Translate generic keywords to their corresponding Orca keywords.
        """

        def inithess():
            """
            Generate an seperate file containing the initial Hessian matrix used as
            guess for the computation.
            """

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
                        ret += ' '.join('{:10.6f}'.format(hess[v + 6 * i][j]) for v in range(n_columns))
                        ret += '\n'
                return ret

            # Check Hessian dimension
            dim = len(value)
            if len(value.shape) == 1:
                dim = int(dim ** 0.5)
                hess = np.reshape(value, (dim, dim))
            else:
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
            hess_path = builtins.config.jm.workdir + "/tmp_hessian.txt"
            with open(hess_path, "w") as hess_file:
                hess_file.write(hess_str)

            settings.specific.orca.geom.InHess = "read"
            settings.specific.orca.geom.InHessName = '"' + hess_path + '"'

            return settings

        def constraint():
            cons = ''
            if isinstance(value, Settings):
                for k, v in value.items():
                    ks = k.split()
                    if ks[0] == 'dist' and len(ks) == 3:
                        cons += '{{ B {:s} {:s} {:f} C }}'.format(*ks[1:], v)
                    elif ks[0] == 'angle' and len(ks) == 4:
                        cons += '{{ A {:s} {:s} {:s} {:f} C }}'.format(*ks[1:], v)
                    elif ks[0] == 'dihed' and len(ks) == 5:
                        cons += '{{ D {:s} {:s} {:s} {:s} {:f} C }}'.format(*ks[1:], v)
                    else:
                        warn('Invalid constraint key: ' + k)
            settings.specific.orca.geom.Constraints._end = cons

        def freeze():
            cons = ''
            for a in value:
                cons += '{{ C {:d} C }}'.format(a)
            settings.specific.orca.geom.Constraints._end = cons

        def selected_atoms():
            settings.specific.dftb.geometry.optim = "cartesian"
            if not isinstance(value, list):
                msg = 'selected_atoms ' + str(value) + ' is not a list'
                raise RuntimeError(msg)
            if isinstance(value[0], int):
                cons = ''
                for a in range(len(mol)):
                    if a not in value:
                        cons += '{{ C {:d} C }}'.format(a)
                settings.specific.orca.geom.Constraints._end = cons
            else:
                cons = ''
                for a in range(len(mol)):
                    if mol.atoms[a].symbol not in value:
                        cons += '{{ C {:d} C }}'.format(a)
                settings.specific.orca.geom.Constraints._end = cons

        # Available translations
        functions = {'inithess': inithess,
                     'freeze': freeze,
                     'selected_atoms': selected_atoms,
                     'constraint': constraint}
        if key in functions:
            functions[key]()
        else:
            msg = 'Keyword ' + key + ' not implemented for package ORCA'
            warn(msg)


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results"""

    def __init__(self, settings, molecule, job_name, plams_dir=None,
                 status='done'):
        properties = package_properties['orca']
        super().__init__(settings, molecule, job_name=job_name,
                         plams_dir=plams_dir, properties=properties,
                         status=status)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, status):
        plams_dir = archive["plams_dir"].path
        return ORCA_Result(settings, molecule, job_name, plams_dir, status)

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
