"""Settings use to describe the input of a Quantum chemistry package."""

from __future__ import annotations

from functools import wraps
from collections.abc import Callable, Mapping
from typing import NoReturn, Any, TypeVar

from scm import plams

__all__ = ['Settings']

_Self = TypeVar("_Self", bound="Settings")


class Settings(plams.core.settings.Settings):
    """A subclass of :class:`plams.Settings<scm.plams.core.settings.Settings>`.

    The difference with respect to plams' Settings are:

    - :code:`settings['a.b']` is equivalent to :code:`settings['a']['b'] = settings.a.b`
    """

    def __getitem__(self, name: str) -> Any:
        """Implement :meth:`self[name]<object.__getitem__>`."""
        return dict.__getitem__(self, name)

    def __setitem__(self, name: str, value: Any) -> None:
        """Implement :meth:`self[name] = value<object.__setitem__>`.

        Like its counterpart in :class:`dict` but passed dictionaries are converted
        into instances of :class:`type(self)<Settings>`.
        """
        if isinstance(value, dict):
            cls = type(self)
            value = cls(value)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name: str) -> Any:
        """Implement :meth:`del self[name]<object.__delitem__>`."""
        dict.__delitem__(self, name)

    def copy(self: _Self) -> _Self:
        """Create a deep(-ish) copy of this instance.

        All nested settings instances embedded within *self* are copied recursively;
        all other objects set without copying.
        """
        cls = type(self)
        ret = cls()
        for name in self:
            if isinstance(self[name], plams.core.settings.Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret

    def __deepcopy__(self: _Self, _: object) -> _Self:
        """Implement :func:`copy.deepcopy(self)<copy.deepcopy>`.

        Serves as an alias for :meth:`Settings.copy`.
        """
        return self.copy()

    def overlay(self: _Self, other: Mapping[str, Any]) -> _Self:
        """Return new instance of :class:`~qmflows.Settings` that is a copy of this instance updated with *other*."""  # noqa: E501
        ret = self.copy()
        ret.update(other)
        return ret


def _unsuported_operation(meth: Callable[..., Any]) -> Callable[..., NoReturn]:
    """Decorate a method such that it raises a :exc:`TypeError`."""
    @wraps(meth)
    def new_meth(self, *args, **kwargs) -> NoReturn:
        raise TypeError(f"{self.__class__.__name__!r} does not "
                        "support item assignment or deletion")
    return new_meth


def _super_method(meth: Callable[..., '_Settings']) -> Callable[..., '_Settings']:
    """Decorate a method such that it's body is executed by the super-class."""
    @wraps(meth)
    def new_meth(self, *args, **kwargs) -> '_Settings':
        cls = type(self)
        func = getattr(super(), meth.__name__)
        ret = func(*args, **kwargs)
        return cls(ret)
    return new_meth


class _Settings(Settings):
    """Immutable-ish counterpart of :class:`Settings`."""

    def __init__(self, *args, **kwargs):
        """Initialize an instance."""
        cls = type(self)
        set_item = super().__setitem__

        dict.__init__(self, *args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                set_item(k, cls(v))
            elif isinstance(v, list):
                set_item(k, [cls(i) if isinstance(i, dict) else i for i in v])

    def __missing__(self, name: str) -> NoReturn:
        """Raise a :exc:`KeyError`."""
        raise KeyError(name)

    __setitem__ = _unsuported_operation(Settings.__setitem__)
    __delitem__ = _unsuported_operation(Settings.__delitem__)
    clear = _unsuported_operation(Settings.clear)
    popitem = _unsuported_operation(Settings.popitem)
    setdefault = _unsuported_operation(Settings.setdefault)
    update = _unsuported_operation(Settings.update)

    overlay = _super_method(Settings.overlay)
    merge = _super_method(Settings.merge)
    copy = _super_method(Settings.copy)
    flatten = _super_method(Settings.flatten)

    __init__.__doc__ = Settings.__init__.__doc__
