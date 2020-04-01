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

ReturnDict = Dict[str, T]  # This should technically be a :class:`typing.TypedDict`
MakeRec = Callable[[T], ReturnDict]


class SerMolecule(Serialiser):
    """Based on the Plams molecule this class encode and decode the information related to the molecule using the JSON format."""  # noqa: E501

    def __init__(self) -> None:
        super().__init__(Molecule)

    def encode(self, obj: Molecule, make_rec: MakeRec) -> ReturnDict:
        return make_rec(obj.as_dict())

    def decode(self, cls: Type[Molecule], data: Mapping) -> Molecule:
        return cls.from_dict(data)


class SerMol(Serialiser):
    """Based on the RDKit molecule this class encodes and decodes the information related to the molecule using a string."""  # noqa: E501

    def __init__(self) -> None:
        super().__init__(Mol)

    def encode(self, obj: Mol, make_rec: MakeRec) -> ReturnDict:
        return make_rec(base64.b64encode(obj.ToBinary()).decode('ascii'))

    def decode(self, cls: Type[Mol], data: str) -> Mol:
        return cls(base64.b64decode(data.encode('ascii')))


class SerSettings(Serialiser):
    """Class to encode and decode the :class:`~qmflows.Settings` class using its internal dictionary structure."""  # noqa: E501

    def __init__(self) -> None:
        super().__init__(Settings)

    def encode(self, obj: Settings, make_rec: MakeRec) -> ReturnDict:
        return make_rec(obj.as_dict())

    def decode(self, cls: Type[Settings], data: Mapping) -> Settings:
        return cls(data)


class SerNDFrame(Serialiser):
    """Class to encode and decode the :class:`pandas.Series` and :class:`pandas.DataFrame` instances."""  # noqa: E501

    def __init__(self, name: Any = NDFrame) -> None:
        super().__init__(name)

    def encode(self, obj: NDFrame, make_rec: MakeRec) -> ReturnDict:
        return make_rec(obj.to_dict())

    def decode(self, cls: Type[NDFrame], data: Mapping) -> NDFrame:
        return cls(data)
