"""A module containing containing :mod:`yaml` loaders with duplicate key checking.

Index
-----
.. currentmodule:: qmflows.yaml_utils
.. autosummary::
    UniqueUnsafeLoader
    UniqueFullLoader
    UniqueSafeLoader

API
---
.. autoclass:: UniqueUnsafeLoader
.. autoclass:: UniqueFullLoader
.. autoclass:: UniqueSafeLoader

"""

from __future__ import annotations

from collections.abc import Hashable
from typing import Any

# Use the fast C-based loaders if possible
try:
    from yaml import (
        CUnsafeLoader as UnsafeLoader,
        CFullLoader as FullLoader,
        CSafeLoader as SafeLoader,
    )
except ImportError:
    from yaml import UnsafeLoader, FullLoader, SafeLoader  # type: ignore[assignment]

from yaml.nodes import MappingNode
from yaml.constructor import ConstructorError, BaseConstructor, SafeConstructor

__all__ = ['UniqueLoader', 'UniqueUnsafeLoader', 'UniqueFullLoader', 'UniqueSafeLoader']


def _construct_mapping(
    loader: BaseConstructor,
    node: MappingNode,
    deep: bool = False,
) -> dict[Any, Any]:
    """A helper function for handling :meth:`~yaml.BaseConstructor.construct_mapping` methods."""
    if not isinstance(node, MappingNode):
        raise ConstructorError(
            None, None,
            f"expected a mapping node, but found {node.id}", node.start_mark,
        )

    if isinstance(loader, SafeConstructor):
        loader.flatten_mapping(node)

    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, Hashable):
            raise ConstructorError(
                "while constructing a mapping", node.start_mark,
                "found unhashable key", key_node.start_mark,
            )
        elif key in mapping:
            raise ConstructorError(
                "while constructing a mapping", node.start_mark,
                "found duplicate key", key_node.start_mark,
            )

        value = loader.construct_object(value_node, deep=deep)
        mapping[key] = value
    return mapping


class UniqueUnsafeLoader(UnsafeLoader):
    """A :class:`~yaml.UnsafeLoader` subclass with duplicate key checking."""

    def construct_mapping(self, node: MappingNode, deep: bool = False) -> dict[Any, Any]:
        """Construct Convert the passed **node** into a :class:`dict`."""
        return _construct_mapping(self, node, deep)


class UniqueFullLoader(FullLoader):
    """A :class:`~yaml.FullLoader` subclass with duplicate key checking."""

    def construct_mapping(self, node: MappingNode, deep: bool = False) -> dict[Any, Any]:
        """Construct Convert the passed **node** into a :class:`dict`."""
        return _construct_mapping(self, node, deep)


class UniqueSafeLoader(SafeLoader):
    """A :class:`~yaml.SafeLoader` subclass with duplicate key checking."""

    def construct_mapping(self, node: MappingNode, deep: bool = False) -> dict[Any, Any]:
        """Construct Convert the passed **node** into a :class:`dict`."""
        return _construct_mapping(self, node, deep)


UniqueLoader = UniqueUnsafeLoader
