"""Settings use to describe the input of a Quantum chemistry package."""

__all__ = ['Settings']

from functools import wraps
from typing import NoReturn, Any, Callable

from scm import plams


class Settings(plams.core.settings.Settings, ):
    """A subclass of :class:`plams.Settings<scm.plams.core.settings.Settings>`.

    The difference with respect to plams' Settings are:
    
    - :code:`settings['a.b']` is equivalent to :code:`settings['a']['b'] = settings.a.b`
    - in :meth:`Settings.update`: :code:`settings.__block_replace = True` results in removal 
      of all existing key value pairs.
      ``__block_replace`` can be either in the updated settings or in the updating settings object.
    """

    def __getitem__(self, name):
        """Implement :meth:`self[name]<object.__getitem__>`."""
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """Implement :meth:`self[name] = value<object.__setitem__>`.

        Like its counterpart in :class:`dict` but passed dictionaries are converted
        into instances of :class:`type(self)<Settings>`.
        """
        if isinstance(value, dict):
            cls = type(self)
            value = cls(value)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """Implement :meth:`del self[name]<object.__delitem__>`."""
        dict.__delitem__(self, name)

    def copy(self):
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

    def __deepcopy__(self, _):
        """Implement :func:`copy.deepcopy(self)<copy.deepcopy>`.

        Serves as an alias for :meth:`Settings.copy`.
        """
        return self.copy()

    def overlay(self, other):
        """Return new instance of |Settings| that is a copy of this instance updated with *other*."""  # noqa: E501
        ret = self.copy()
        ret.update(other)
        return ret

    def update(self, other):
        """Implement :meth:`self.update(other)<dict.update>`.

        Like PLAMS update, but:
        :code:`__block_replace = True` results in removing all existing key value pairs
        in the current node of the settings tree

        For example:

        .. code:: python

            >>> from qmflows import Settings

            >>> t=Settings()
            >>> t.input.xc.lda = ""

            >>> u = Settings()
            >>> u.input.xc.hybrid = 'b3lyp'
            >>> print(t.overlay(u))
            input:
                xc:
                    hybrid:     b3lyp
                    lda:

            >>> u.input.xc.__block_replace = True
            >>> print(t.overlay(u))
            input:
                xc:
                    __block_replace:     True
                    hybrid:     b3lyp

        """
        cls = type(self)
        for name in other:
            if isinstance(other[name], cls):
                if name not in self or not isinstance(self[name], cls):
                    self[name] = other[name].copy()
                else:
                    # _block_replace can be used to remove all existing key value pairs
                    # in an input block
                    br = '__block_replace'
                    pred1 = (br in other[name] and other[name][br])
                    pred2 = (br in self[name] and self[name][br])
                    if pred1 or pred2:
                        self[name] = cls()
                    self[name].update(other[name])
            else:
                self[name] = other[name]


def _unsuported_operation(meth: Callable) -> Callable[..., NoReturn]:
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
