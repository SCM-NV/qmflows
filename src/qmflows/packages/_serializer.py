"""Various serialisers used in QMFlows."""

from __future__ import annotations

import base64
from collections.abc import Callable, Iterable, Sequence
from typing import TypeVar, Any, TYPE_CHECKING

from noodles.serial import Serialiser

if TYPE_CHECKING:
    from .._settings import Settings
    from scm.plams import Molecule
    from rdkit.Chem import Mol
    from pandas.core.generic import NDFrame
    from typing import Protocol, TypedDict, Generic
    from typing_extensions import NotRequired

    T = TypeVar('T')

    _RecDictBase = TypedDict('_RecDictBase', {'type': str, 'class': str})

    class _RecDict(_RecDictBase, Generic[T]):
        _noodles: str
        data: T
        ref: NotRequired[bool]
        host: NotRequired[None | str]
        files: NotRequired[list[str]]

    class _MakeRec(Protocol[T]):
        def __call__(
                self,
                __data: T,
                *,
                ref: bool = ...,
                files: list[str] = ...,
            ) -> _RecDict[T]: ...

    class _AtomDict(TypedDict):
        atnum: int
        bonds: list[int]
        coords: tuple[float, float, float]
        properties: Settings

    class _BondDict(TypedDict):
        atom1: int
        atom2: int
        order: float
        properties: Settings

    class _MolDict(TypedDict):
        charge: NotRequired[int]
        atoms: list[_AtomDict]
        bonds: list[_BondDict]
        lattice: list[Sequence[float]]
        properties: Settings

__all__ = ['SerMolecule', 'SerMol', 'SerSettings', 'SerNDFrame', 'SerReduce']


class SerMolecule(Serialiser):
    """Based on the Plams molecule this class encode and decode the information related to the molecule using the JSON format."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerMolecule` instance."""
        super().__init__("Molecule")

    def encode(self, obj: Molecule, make_rec: _MakeRec[_MolDict]) -> _RecDict[_MolDict]:
        """Encode the passed PLAMS Molecule."""
        return make_rec(obj.as_dict())

    def decode(self, cls: type[Molecule], data: _MolDict) -> Molecule:
        """Decode the passed data into a PLAMS Molecule."""
        return cls.from_dict(data)


class SerMol(Serialiser):
    """Based on the RDKit molecule this class encodes and decodes the information related to the molecule using a string."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerMol` instance."""
        super().__init__("Mol")

    def encode(self, obj: Mol, make_rec: _MakeRec[str]) -> _RecDict[str]:
        """Encode the passed RDKit Mol."""
        mol_bytes: bytes = obj.ToBinary()
        return make_rec(base64.b64encode(mol_bytes).decode('ascii'))

    def decode(self, cls: type[Mol], data: str) -> Mol:
        """Decode the passed data into a RDKit Mol."""
        return cls(base64.b64decode(data.encode('ascii')))


class SerSettings(Serialiser):
    """Class to encode and decode the :class:`~qmflows.Settings` class using its internal dictionary structure."""  # noqa: E501

    def __init__(self) -> None:
        """Initialize a :class:`SerSettings` instance."""
        super().__init__("Settings")

    def encode(self, obj: Settings, make_rec: _MakeRec[dict[Any, Any]]) -> _RecDict[dict[Any, Any]]:
        """Encode the passed PLAMS Settings."""
        return make_rec(obj.as_dict())

    def decode(self, cls: type[Settings], data: dict[Any, Any]) -> Settings:
        """Decode the passed data into a PLAMS Settings."""
        return cls(data)


class SerNDFrame(Serialiser):
    """Class to encode and decode the :class:`pandas.Series` and :class:`pandas.DataFrame` instances."""  # noqa: E501

    def __init__(self, name: Any = "NDFrame") -> None:
        """Initialize a :class:`SerNDFrame` instance."""
        super().__init__(name)

    def encode(self, obj: NDFrame, make_rec: _MakeRec[dict[str, Any]]) -> _RecDict[dict[str, Any]]:
        """Encode the passed pandas Series or DataFrame."""
        return make_rec(obj.to_dict())

    def decode(self, cls: type[NDFrame], data: dict[str, Any]) -> NDFrame:
        """Decode the passed data into a pandas Series or DataFrame."""
        return cls(data)


class SerReduce(Serialiser):
    """Class to encode :meth:`object.__reduce__`-able objects."""

    def __init__(self, name: Any) -> None:
        """Initialize a :class:`SerReduce` instance."""
        super().__init__(name)

    def encode(
        self, obj: Any, make_rec: Callable[[Any], Any]
    ) -> tuple[tuple[Any, ...], None | Any]:
        """Encode the passed reduce-able object."""
        _, args, *tail = obj.__reduce__()
        state = tail[0] if len(tail) else None
        return args, state

    def decode(self, cls: type[T], data: tuple[Iterable[Any], None | Any]) -> T:
        """Decode the passed data into a PLAMS Molecule."""
        args, state = data
        ret = cls(*args)
        if state is not None and hasattr(ret, "__setstate__"):
            ret.__setstate__(state)
        return ret
