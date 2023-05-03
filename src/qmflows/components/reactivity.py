"""A module containing the :class:`Coordinate` baseclass and its subclasses."""

from __future__ import annotations

from collections.abc import Callable
from typing import ClassVar, overload

import scm.plams.interfaces.molecule.rdkit as molkit
from scm.plams import Molecule
from rdkit.Chem import AllChem

from ..type_hints import MolType
from .._settings import Settings

__all__ = ['Distance', 'Angle', 'Dihedral']


class Coordinate:
    """The :class:`Coordinate` base class."""

    fmt: ClassVar[str] = "{}"
    atoms: tuple[int, ...]

    @property
    def fun(self) -> Callable[..., float]:
        """Getter and setter for the :attr:`Coordinate.fun` property.

        Setting will simply assign the value.
        Getting will return the value and, if it has not been set,
        raise a :exc:`NotImplementedError`.

        """
        try:
            return self._fun
        except AttributeError as ex:
            msg = f'method {self.__class__.__name__}.fun() is not implemented'
            raise NotImplementedError(msg) from ex

    @fun.setter
    def fun(self, value: Callable[..., float]) -> None:
        self._fun = value

    def __init__(self, *args: int) -> None:
        """Initialize a :class:`Coordinate` instance."""
        self.atoms = args

    def get_current_value(self, mol: MolType) -> float:
        """Return the value of the coordinate."""
        if isinstance(mol, Molecule):
            mol = molkit.to_rdmol(mol)
        conf = mol.GetConformer()

        # list of indices
        xs = [i - 1 for i in self.atoms]
        return self.fun(conf, *xs)

    @overload
    def get_settings(self, value: float, mol: None) -> Settings: ...
    @overload
    def get_settings(self, value: None, mol: MolType) -> Settings: ...

    def get_settings(self, value=None, mol=None):
        """Map a :class:`str` representation of :attr:`Coordinate.atoms` to *value*."""
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
    """Class defining an atomic distance."""

    fun: Callable[[AllChem.Conformer, int, int], float]
    fmt: ClassVar[str] = "dist {:d} {:d}"

    def __init__(self, atom1: int, atom2: int) -> None:
        """Initialize a :class:`Distance` instance."""
        super().__init__(atom1, atom2)
        self.fun = AllChem.GetBondLength


class Angle(Coordinate):
    """Class defining an atomic angle."""

    fun: Callable[[AllChem.Conformer, int, int, int], float]
    fmt: ClassVar[str] = "angle {:d} {:d} {:d}"

    def __init__(self, atom1: int, atom2: int, atom3: int) -> None:
        """Initialize an :class:`Angle` instance."""
        super().__init__(atom1, atom2, atom3)

    def get_current_value(self, mol: MolType,
                          rad: bool = False) -> float:
        """Return the value of the coordinate."""
        self.fun = AllChem.GetAngleRad if rad else AllChem.GetAngleDeg
        return super().get_current_value(mol)


class Dihedral(Coordinate):
    """Class defining an atomic dihedral angle."""

    fun: Callable[[AllChem.Conformer, int, int, int, int], float]
    fmt: ClassVar[str] = "dihed {:d} {:d} {:d} {:d}"

    def __init__(self, atom1: int, atom2: int, atom3: int, atom4: int) -> None:
        """Initialize a :class:`Dihedral` instance."""
        super().__init__(atom1, atom2, atom3, atom4)

    def get_current_value(self, mol: MolType,
                          rad: bool = False) -> float:
        """Return the value of the coordinate."""
        self.fun = AllChem.GetDihedralRad if rad else AllChem.GetDihedralDeg
        return super().get_current_value(mol)
