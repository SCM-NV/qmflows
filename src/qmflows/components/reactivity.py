
__all__ = ['Distance', 'Angle', 'Dihedral', 'PES']

from qmflows import templates
from qmflows.settings import Settings
from noodles import gather, schedule
from scm.plams import Molecule
from rdkit.Chem import AllChem
import scm.plams.interfaces.molecule.rdkit as molkit


class Coordinate:
    def __init__(self, *args):
        self.atoms = args
        self.fmt = "{}"
        self.fun = None

    def get_current_value(self, mol):
        """
        Value of the coordinate
        """
        if isinstance(mol, Molecule):
            mol = molkit.to_rdmol(mol)
        conf = mol.GetConformer()

        # list of indices
        xs = [i - 1 for i in self.atoms]
        return self.fun(conf, *xs)

    def get_settings(self, value=None, mol=None):
        s = Settings()
        if value is None and mol is None:
            msg = 'coordinate constraint settings requires a value or molecule'
            raise RuntimeError(msg)
        elif value is None:
            value = self.get_current_value(mol)

        # create settings entry
        data = self.fmt.format(*self.atoms)
        s[data] = value
        return s


class Distance(Coordinate):
    """
    Class defining an atomic distance
    """
    def __init__(self, atom1, atom2):
        super().__init__(atom1, atom2)
        self.fmt = "dist {:d} {:d}"
        self.fun = AllChem.GetBondLength


class Angle(Coordinate):
    """
    Class defining an atomic angle
    """
    def __init__(self, atom1, atom2, atom3):
        super().__init__(atom1, atom2, atom3)
        self.fmt = "angle {:d} {:d} {:d}"

    def get_current_value(self, mol, rad=False):
        self.fun = AllChem.GetAngleRad if rad else AllChem.GetAngleDeg
        return super().get_current_value(mol)


class Dihedral(Coordinate):
    """
    Class defining an atomic dihedral angle
    """
    def __init__(self, atom1, atom2, atom3, atom4):
        super().__init__(atom1, atom2, atom3, atom4)
        self.fmt = "dihed {:d} {:d} {:d} {:d}"

    def get_current_value(self, mol, rad=False):
        self.fun = AllChem.GetDihedralRad if rad else AllChem.GetDihedralDeg
        return super().get_current_value(mol)


@schedule
class PES:
    def __init__(
            self, molecule=None, constraints=None, offset=None,
            get_current_values=False,
            nsteps=0, stepsize=0.0, nested_PES=None):

        self.molecule = molkit.to_rdmol(molecule)
        self.constraints = constraints
        if isinstance(constraints, list):
            self.start = []
            for i, cs in enumerate(constraints):
                if offset is None:
                    self.start.append(0.0)
                else:
                    self.start.append(offset[i])
                if get_current_values:
                    self.start[i] += cs.get_current_value(self.molecule)
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
