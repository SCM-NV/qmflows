"""A module with type-hint related objects used throughout QMFlows.

Index
-----
.. currentmodule:: qmflows.type_hints
.. autosummary::
    WarnMap
    WarnDict
    WarnParser
    MappingScalar
    MappingSequence
    Generic2Special
    PathLike

API
---
.. autoclass:: ParseWarning

.. autodata:: WarnMap
    :annotation: = Mapping[str, ParseWarning]

.. autodata:: WarnDict
    :annotation: = Dict[str, Type[Warning]]

.. autodata:: WarnParser
    :annotation: = Callable[[str, WarnMap], Optional[WarnDict]]

.. autodata:: WarnParseMap
    :annotation: = Mapping[Type[Warning], Callable[[str], Optional[str]]]

.. autodata:: MappingScalar
    :annotation: = MutableMapping[str, Union[None, str, float]]

.. autodata:: MappingSequence
    :annotation: = MutableMapping[str, Union[Sequence[Optional[str]], Sequence[float]]]

.. autodata:: Generic2Special
    :annotation: = Callable[[qmflows.Settings, str, Any, plams.Molecule], None]

.. autodata:: PathLike
    :annotation: = Union[str, bytes, os.PathLike]

"""

from typing import (Mapping, Dict, Type, Callable, Optional, MutableMapping, Union,
                    Sequence, Any, AnyStr, TYPE_CHECKING)

from .common import ParseWarning

if TYPE_CHECKING:
    from os import PathLike as _PathLike
    from scm.plams import Molecule
    from .settings import Settings
else:
    _PathLike = 'os.PathLike'
    Molecule = 'scm.plams.mol.molecule.Molecule'
    Settings = 'qmflows.settings.Settings'

__all__ = [
    'WarnMap', 'WarnDict', 'WarnParser',
    'MappingScalar', 'MappingSequence',
    'Generic2Special',
    'PathLike'
]


#: A Mapping which maps a :class:`str` to a :class:`ParseWarning` instance.
WarnMap = Mapping[str, ParseWarning]

#: A dictionary which maps a :class:`str` to a :exc:`Warning` type.
WarnDict = Dict[str, Type[Warning]]

#: A callable which takes a :class:`str` and :data:`WarnMap` as argument and
#: returns either ``None`` or a :data:`WarnDict`.
WarnParser = Callable[[str, WarnMap], Optional[WarnDict]]

#: A mutable mapping which maps a :class:`str` to a scalar.
MappingScalar = MutableMapping[str, Union[None, str, float]]

#: A mutable mapping which maps a :class:`str` to
#: a :class:`~collections.abc.Sequence` of scalars.
MappingSequence = MutableMapping[str, Union[Sequence[Optional[str]], Sequence[float]]]

#: A function for converting generic to specific settings.
Generic2Special = Callable[[Settings, str, Any, Molecule], None]

#: An alias for all path-like objects.
#: Includes the likes of :class:`str`, :class:`bytes` and :class:`pathlib.Path`.
PathLike = Union[AnyStr, _PathLike]
