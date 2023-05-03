"""Check yaml to Settings conversion."""

from __future__ import annotations

from collections.abc import Callable
from typing import overload, Any

import yaml

from ._settings import Settings
from .type_hints import T
from .yaml_utils import UniqueSafeLoader

__all__ = ['yaml2Settings']


@overload
def yaml2Settings(xs: str | bytes) -> Settings: ...
@overload
def yaml2Settings(xs: str | bytes, mapping_type: Callable[[dict[str, Any]], T]) -> T: ...


def yaml2Settings(xs, mapping_type=Settings):
    """Transform a string containing some data in .yaml format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()

    dct = yaml.load(xs, Loader=UniqueSafeLoader)  # yaml object must be string
    return mapping_type(dct)
