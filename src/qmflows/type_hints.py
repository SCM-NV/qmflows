"""A module with type-hint related objects used throughout QMFlows."""

from typing import Mapping, Dict, Type, Callable, Optional, MutableMapping, Union, Sequence, Any

from scm import plams

from .settings import Settings

__all__ = ['WarnMap', 'WarnDict', 'WarnParser', 'MappingScalar', 'MappingSequence',
           'Generic2Special']


#: A Mapping mapping a str to a Warning type.
WarnMap = Mapping[str, Type[Warning]]

#: A dictionary mapping a str to a Warning type.
WarnDict = Dict[str, Type[Warning]]

#: A callable which takes a str and WarnMap as argument and
#: returns either ``None`` or a WarnDict.
WarnParser = Callable[[str, WarnMap], Optional[WarnDict]]

#: A mutable mapping which maps a str to a scalar.
MappingScalar = MutableMapping[str, Union[None, str, float]]

#: A mutable mapping which maps a str to a sequence of scalars.
MappingSequence = MutableMapping[str, Sequence[Union[None, str, float]]]

#: A function for converting generic to specific settings.
Generic2Special = Callable[[Settings, str, Any, plams.Molecule], None]
