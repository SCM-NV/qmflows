"""Settings use to describe the input of a Quantum chemistry package."""

__all__ = ['Settings']

from functools import wraps
from typing import NoReturn, Any, Callable

from scm import plams


class Settings(plams.core.settings.Settings, ):
    """This is a subclass of the :class:`plams.core.settings.Settings`.

    The difference with respect to plams' Settings are:
    - settings['a.b'] is equivalent to settings['a']['b'] = settings.a.b
    - in update(): settings.__block_replace = True results in removal of all
    existing key value pairs.

      __block_replace can be either in the updated settings or in the updating settings object.
    """

    def __getitem__(self, name):
        """
        Like plams Settings.__getitem__, but
        "settings['a.b'] == 'c'" is equivalent to "settings['a']['b'] == 'c'"
        """
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """
        Like plams Settings.__setitem__, but
        "settings['a.b'] = 'c'" is equivalent to "settings['a']['b'] = 'c'"
        """
        cls = type(self)
        if isinstance(value, dict):
            value = cls(value)
        dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """
        Like plams Settings.__setitem__, but
        "del settings['a.b']" is equivalent to "del settings['a']['b'] = 'c'"
        """
        dict.__delitem__(self, name)

    def copy(self):
        cls = type(self)
        ret = cls()
        for name in self:
            if isinstance(self[name], plams.core.settings.Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret

    def __deepcopy__(self, _):
        return self.copy()

    def overlay(self, other):
        """
        Return new instance of |Settings| that is a copy of this instance
        updated with *other*.
        """
        ret = self.copy()
        ret.update(other)
        return ret

    def update(self, other):
        """
        Like PLAMS update, but:
        __block_replace = True results in removing all existing key value pairs
        in the current node of the settings tree

        For example:

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


class _Settings(Settings):
    """Immutable-ish counterpart of :class:`Settings`."""

    def __init__(self, *args, **kwargs):
        cls = type(self)
        set_item = super().__setitem__

        dict.__init__(self, *args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                set_item(k, cls(v))
            elif isinstance(v, list):
                set_item(k, [cls(i) if isinstance(i, dict) else i for i in v])

    @staticmethod
    def unsuported_operation(meth: Callable) -> Callable[..., NoReturn]:
        """Decorate a method such that it raises a :exc:`TypeError`."""
        @wraps(meth)
        def new_meth(self, *args, **kwargs) -> NoReturn:
            raise TypeError(f"{self.__class__.__name__!r} does not "
                            "support item assignment or deletion")
        return new_meth

    @staticmethod
    def super_method(meth: Callable[..., '_Settings']) -> Callable[..., '_Settings']:
        """Decorate a method such that it's body is executed by the super-class."""
        @wraps(meth)
        def new_meth(self, *args, **kwargs) -> '_Settings':
            cls = type(self)
            func = getattr(super(), meth.__name__)
            ret = func(*args, **kwargs)
            return cls(ret)
        return new_meth

    __setitem__ = unsuported_operation(Settings.__setitem__)
    __delitem__ = unsuported_operation(Settings.__delitem__)
    clear = unsuported_operation(Settings.clear)
    popitem = unsuported_operation(Settings.popitem)
    setdefault = unsuported_operation(Settings.setdefault)
    update = unsuported_operation(Settings.update)

    overlay = super_method(Settings.overlay)
    merge = super_method(Settings.merge)
    copy = super_method(Settings.copy)
    flatten = super_method(Settings.flatten)
