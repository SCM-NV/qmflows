"""Check yaml to Settings conversion."""

__all__ = ['yaml2Settings']

import yaml
from typing import AnyStr
from .settings import Settings


def yaml2Settings(xs: AnyStr) -> Settings:
    """Transform a string containing some data in .yaml format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = yaml.load(xs, Loader=yaml.FullLoader)  # yaml object must be string
    return Settings(s)
