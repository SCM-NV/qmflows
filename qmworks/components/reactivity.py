
__all__ = ['Distance', 'Angle', 'Dihedral', 'PES']

from qmworks import templates, molkit
from qmworks.settings import Settings
from noodles import gather, schedule
from plams import Molecule
from rdkit.Chem import AllChem


class Distance:
    """
    Class defining an atomic distance
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def get_current_value(self, mol):
        if isinstance(mol, Molecule):
            mol = molkit.to_rdmol(mol)
        conf = mol.GetConformer()
        return AllChem.GetBondLength(conf, self.atom1, self.atom2)

    def get_settings(self, value=None, mol=None):
        s = Settings()
        if value is None:
            if mol is None:
                msg = 'Distance constraint settings requires a value or molecule'
                raise RuntimeError(msg)
            else:
                value = self.get_current_value(mol)
        s["dist {:d} {:d}".format(self.atom1, self.atom2)] = value
        return s


class Angle:
    """
    Class defining an atomic angle
    """
    def __init__(self, atom1, atom2, atom3):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def get_current_value(self, mol, rad=False):
        if isinstance(mol, Molecule):
            mol = molkit.to_rdmol(mol)
        conf = mol.GetConformer()
        if rad:
            return AllChem.GetAngleRad(conf, self.atom1, self.atom2, self.atom3)
        else:
            return AllChem.GetAngleDeg(conf, self.atom1, self.atom2, self.atom3)

    def get_settings(self, value=None, mol=None):
        s = Settings()
        if value is None:
            if mol is None:
                msg = 'Angle constraint settings requires a value or molecule'
                raise RuntimeError(msg)
            else:
                value = self.get_current_value(mol)
        s["angle {:d} {:d} {:d}".format(self.atom1, self.atom2, self.atom3)] = value
        return s


class Dihedral:
    """
    Class defining an atomic dihedral angle
    """
    def __init__(self, atom1, atom2, atom3, atom4):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def get_current_value(self, mol, rad=False):
        if isinstance(mol, Molecule):
            mol = molkit.to_rdmol(mol)
        conf = mol.GetConformer()
        if rad:
            return AllChem.GetDihedralRad(conf, self.atom1, self.atom2, self.atom3, self.atom4)
        else:
            return AllChem.GetDihedralDeg(conf, self.atom1, self.atom2, self.atom3, self.atom4)

    def get_settings(self, value=None, mol=None):
        s = Settings()
        if value is None:
            if mol is None:
                msg = 'Dihedral constraint settings requires a value or molecule'
                raise RuntimeError(msg)
            else:
                value = self.get_current_value(mol)
        s["dihed {:d} {:d} {:d} {:d}".format(self.atom1, self.atom2, self.atom3, self.atom4)] = value
        return s


@schedule
class PES:
    def __init__(self, molecule=None, constraints=None, offset=None, get_current_values=False,
                 nsteps=0, stepsize=0.0, nested_PES=None):
        self.molecule = molkit.to_rdmol(molecule)
        self.constraints = constraints
        if isinstance(constraints, list):
            self.start = []
            for i in range(len(constraints)):
                if offset is None:
                    self.start.append(0.0)
                else:
                    self.start.append(offset[i])
                if get_current_values:
                    self.start[i] += constraints[i].get_current_value(self.molecule)
        else:
            if offset is None:
                self.start = 0.0
            else:
                self.start = offset
            if get_current_values:
                self.start += constraints.get_current_value(self.molecule)
        self.nsteps = nsteps
        self.stepsize = stepsize
        self.nested_PES = nested_PES

    def scan(self, package, settings, job_name="PESscan"):
        """
        This function enables multidimensional potential energy surface scans.
        The argument 'scan' should be a dictionary with keys 'constraint' and 'surface'
        and optionally 'scan'. The latter allows nested (multidimensional) scans.

        Example
        .. code-block:: python

            scan = {'constraint': "dist 1 5",
                    'surface': {'nsteps':6, 'start': 2.3, 'stepsize': 0.1},
                    'scan': {'constraint': "dist 3 4",
                             'surface': {'nsteps':6, 'start': 2.3, 'stepsize': 0.1}
                             }
                   }

        It is also possible to scan over two coordinates in a concerted way, for example:
        .. code-block:: python

            scan = {'constraint': ["dist 1 5", "dist 3 4"],
                       'surface': {'nsteps':6, 'start':[2.3,2.3], 'stepsize': [0.1,0.1]} }

        In this example multiple listed values for 'start' and 'stepsize' correspond to
        the same number of listed 'constraint's

        :returns:
            A list of promised result objects
        """
        pes_jobs = []

        for step in range(self.nsteps):
            new_settings = settings.overlay(self.get_constraint_settings(step))
            name = job_name + "_" + str(step)
            if self.nested_PES:
                pes_jobs.append(self.nested_PES.scan(package, new_settings, name))
            else:
                pes_jobs.append(self.pes_job(package, new_settings, name))
        return gather(*pes_jobs)

    def pes_job(self, package, settings, job_name):
        if isinstance(package, list):
            name = job_name + "_opt"
            optimized_mol = package[0](templates.geometry.overlay(settings),
                                       self.molecule, job_name=name).molecule
            result = package[1](templates.singlepoint.overlay(settings),
                                optimized_mol, job_name=job_name)
        else:
            result = package(templates.geometry.overlay(settings), self.molecule,
                             job_name=job_name)
        return result

    def get_constraint_settings(self, step):
        s = Settings()
        if isinstance(self.constraints, list):
            for c in range(len(self.constraints)):
                s.constraint.update(
                    self.constraints[c].get_settings(self.start[c] +
                                                     self.stepsize[c] * step))
        else:
            s.constraint = self.constraints.get_settings(self.start + self.stepsize * step)
        return s
