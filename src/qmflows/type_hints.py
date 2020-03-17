"""A module with type-hint related objects used throughout QMFlows.

Index
-----
.. currentmodule:: qmflows.type_hints
.. autosummary::
    {autosummary}

API
---
{autodata}

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


__doc__ = __doc__.format(autosummary='\n'.join(f'    {i}' for i in __all__),
                         autodata='\n'.join(f'.. autodata:: {i}' for i in __all__))
