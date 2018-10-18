
__all__ = ['Distance', 'Angle', 'Dihedral']

from qmflows.settings import Settings
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
