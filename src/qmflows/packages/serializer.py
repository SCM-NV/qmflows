"""Various serialisers used in QMFlows."""

import base64
from typing import Type, TypeVar, Callable, Any, Dict, Mapping, TYPE_CHECKING

from noodles.serial import Serialiser

if TYPE_CHECKING:
    from ..settings import Settings
    from scm.plams import Molecule
    from rdkit.Chem import Mol
    from pandas.core.generic import NDFrame
    from pandas import DataFrame, Series

else:  # Don't bother importing all this stuff when not type checking
    Settings = 'qmflows.settings.Settings'
    Molecule = 'scm.plams.mol.molecule.Molecule'
    Mol = 'rdkit.Chem.rdchem.Mol'
    NDFrame = 'pandas.core.generic.NDFrame'
    DataFrame = 'pandas.core.frame.DataFrame'
    Series = 'pandas.core.series.Series'

__all__ = ['SerMolecule', 'SerMol', 'SerSettings', 'SerNDFrame']

T = TypeVar('T')


class SerMolecule(Serialiser):
    """Based on the Plams molecule this class encode and decode the information related to the molecule using the JSON format."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerMolecule` instance."""
        super().__init__(Molecule)

    def encode(self, obj: Molecule, make_rec: Callable[[T], Dict[str, T]]) -> Dict[str, T]:
        """Encode the passed PLAMS Molecule."""
        return make_rec(obj.as_dict())

    def decode(self, cls: Type[Molecule], data: Mapping) -> Molecule:
        """Decode the passed data into a PLAMS Molecule."""
        return cls.from_dict(data)


class SerMol(Serialiser):
    """Based on the RDKit molecule this class encodes and decodes the information related to the molecule using a string."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerMol` instance."""
        super().__init__(Mol)

    def encode(self, obj: Mol, make_rec: Callable[[T], Dict[str, T]]) -> Dict[str, T]:
        """Encode the passed RDKit Mol."""
        return make_rec(base64.b64encode(obj.ToBinary()).decode('ascii'))

    def decode(self, cls: Type[Mol], data: str) -> Mol:
        """Decode the passed data into a RDKit Mol."""
        return cls(base64.b64decode(data.encode('ascii')))


class SerSettings(Serialiser):
    """Class to encode and decode the :class:`~qmflows.Settings` class using its internal dictionary structure."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerSettings` instance."""
        super().__init__(Settings)

    def encode(self, obj: Settings, make_rec: Callable[[T], Dict[str, T]]) -> Dict[str, T]:
        """Encode the passed PLAMS Settings."""
        return make_rec(obj.as_dict())

    def decode(self, cls: Type[Settings], data: Mapping) -> Settings:
        """Decode the passed data into a PLAMS Settings."""
        return cls(data)


class SerNDFrame(Serialiser):
    """Class to encode and decode the :class:`pandas.Series` and :class:`pandas.DataFrame` instances."""  # noqa: E501

    def __init__(self, name: Any = NDFrame) -> None:
        """Initialize a :class:`SerNDFrame` instance."""
        super().__init__(name)

    def encode(self, obj: NDFrame, make_rec: Callable[[T], Dict[str, T]]) -> Dict[str, T]:
        """Encode the passed pandas Series or DataFrame."""
        return make_rec(obj.to_dict())

    def decode(self, cls: Type[NDFrame], data: Mapping) -> NDFrame:
        """Decode the passed data into a pandas Series or DataFrame."""
        return cls(data)
