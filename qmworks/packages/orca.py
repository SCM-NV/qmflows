# =======>  Standard and third party Python Libraries <======
from noodles import (Storable)
from os.path import join
from qmworks.settings import Settings
from qmworks.packages.packages import (Package, Result)
from qmworks.parsers.orca_parser import parse_molecule

import builtins
import plams
import numpy as np
# ============================= Orca ==========================================


class _PropertyGetter(Storable):
    def __init__(self, result, section):
        super(_PropertyGetter, self).__init__()
        self.result = result
        self.section = section

    def __getattr__(self, prop):
        pass

    def __getitem__(self, prop):
        pass


class ORCA(Package):
    """
    This class prepare the input to run a Orca job using both Plams and
    templates. It also does the manangement of the input/output files resulting
    from running Orca and returns a Results object that containing the methods
    and data required to retrieve the output.
    """
    def __init__(self):
        super(ORCA, self).__init__("orca")
        self.generic_dict_file = 'generic2ORCA.json'

    def prerun(self):
        pass

    def run_job(self, settings, mol, job_name="ORCAjob"):

        orca_settings = Settings()
        orca_settings.input = settings.specific.orca

        result = plams.ORCAJob(molecule=mol, settings=orca_settings,
                               name=job_name).run()

        return ORCA_Result(orca_settings, mol, result.job.name, result.job.path)

    def postrun(self):
        pass

    def handle_special_keywords(self, settings, key, value, mol):
        if key == "inithess":
            hess = value
            d = len(value)
            if len(value.shape) == 1:
                d = len(value)**0.5
                if int(d) == d:
                    d = int(d)
                    hess = np.reshape(value, (d, d))
            hess_path = builtins.config.jm.workdir + "/tmp_hessian.txt"
            hess_file = open(hess_path, "w")
            hess_file.write('\n$orca_hessian_file\n\n$hessian\n' + str(d) + '\n')
            for i in range(int((d-1) / 6) + 1):
                hess_file.write("         " + " ".join(['{:10d}'.format(v + 6 * i) for v in range(min(6, d - 6*i))]) + '\n')
                for j in range(len(hess[i])):
                    hess_file.write('{:7d}'.format(j) + "     " + " ".join(['{:10.6f}'.format(hess[6*i+v][j]) for v in range(min(6, d - 6*i))]) + '\n')
            hess_file.write('\n$atoms\n')
            hess_file.write(str(len(mol)) + '\n')
            for a in mol.atoms:
                hess_file.write('{:2s}{:12.4f}{:14.6f}{:14.6f}{:14.6f}\n'.format(a.symbol, a._getmass(), *a.coords))
            hess_file.write('\n\n$end\n')
            hess_file.close()
            settings.specific.orca.geom.InHess = "read"
            settings.specific.orca.geom.InHessName = '"' + hess_path + '"'
        else:
            warn(UserWarning('Keyword ' + key + ' doesn\'t exist'))


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
