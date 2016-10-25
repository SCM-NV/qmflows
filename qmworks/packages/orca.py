# =======>  Standard and third party Python Libraries <======
from os.path import join
from qmworks.settings import Settings
from qmworks.packages.packages import (Package, Result)
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

        result = plams.ORCAJob(molecule=mol, settings=orca_settings,
                               name=job_name).run()

        return ORCA_Result(orca_settings, mol, result.job.name, result.job.path)

    def postrun(self):
        pass

    @staticmethod
    def handle_special_keywords(settings, key, value, mol):
        """
        Translate generic keywords to their corresponding Orca keywords.
        """
        # Available translations
        funs = {'inithess': translate_inithess}
        f = funs.get(key)
        # Apply translation
        if f is not None:
            return f(settings, key, value, mol)
        else:
            msg = 'Keyword ' + key + ' doesn\'t exist'
            warn(msg)


def translate_inithess(settings, key, value, mol):
    """
    Generate an seperate file containing the initial Hessian matrix used as
    guess for the computation.
    """
    def format_atom(at):
        s, m, xs = at.symbol, at._getmass(), at.coords
        return '{:2s}{:12.4f}{:14.6f}{:14.6f}{:14.6f}\n'.format(s, m, *xs)

    def format_hessian(d, hess):
        """ Format numpy array to Orca matrix format """
        acc = ''
        for i in range((d - 1) // 6 + 1):
            m = min(6, d - 6 * i)
            just = 13 + m * 11  # Right justification
            xs = ''.join('{:^11d}'.format(v + 6 * i) for v in range(m))
            acc += xs.rjust(just)  + '\n'
            for j in range(len(hess[i])):
                ys = '{:7d}'.format(j).ljust(11)
                ys += ''.join('{:11.7f}'.format(hess[v + 6 * i][j]) for v in range(m))
                acc += ys + '\n'
        return acc

    # Check Hessian dimension
    d = len(value)
    if len(value.shape) == 1:
        dim = int(d ** 0.5)
        hess = np.reshape(value, (dim, dim))
    else:
        hess = value

    # Header
    hess_str = '\n$orca_hessian_file\n\n$hessian\n' + str(d) + '\n'
    # Actual Hessian
    hess_str += format_hessian(d, hess)
    # Atoms header
    hess_str += '\n$atoms\n'
    # Atoms coordinates
    hess_str += str(len(mol)) + '\n'
    hess_str += ''.join(format_atom(at) for at in mol.atoms)
    # The end
    hess_str += '\n\n$end\n'

    # Store the hessian in the plams_dir
    hess_path = builtins.config.jm.workdir + "/tmp_hessian.txt"
    with open(hess_path, "w") as  hess_file:
        hess_file.write(hess_str)

    settings.specific.orca.geom.InHess = "read"
    settings.specific.orca.geom.InHessName = '"' + hess_path + '"'

    return settings


class ORCA_Result(Result):
    """Class providing access to PLAMS OrcaJob results"""

    def __init__(self, settings, molecule, job_name, plams_dir, project_name=None):
        properties = 'data/dictionaries/propertiesORCA.json'
        super().__init__(settings, molecule, job_name=job_name, plams_dir=plams_dir,
                         project_name=project_name, properties=properties)

    @classmethod
    def from_dict(cls, settings, molecule, job_name, archive, project_name):
        plams_dir = archive["plams_dir"].path
        return ORCA_Result(settings, molecule, job_name, plams_dir, project_name)

    @property
    def molecule(self):
        """ Retrieve the molecule from the output file"""
        plams_dir = self.archive["plams_dir"].path
        file_name = join(plams_dir, '{}.out'.format(self.job_name))
        return parse_molecule(file_name, self._molecule)

orca = ORCA()
