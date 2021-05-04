"""Check yaml to Settings conversion."""

__all__ = ['yaml2Settings']

from typing import Union, Callable

import yaml

from .settings import Settings
from .type_hints import T
from .yaml_utils import UniqueSafeLoader


def yaml2Settings(
    xs: Union[str, bytes],
    mapping_type: Callable[[dict], T] = Settings,
) -> T:
    """Transform a string containing some data in .yaml format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()

    dct = yaml.load(xs, Loader=UniqueSafeLoader)  # yaml object must be string
    return mapping_type(dct)
