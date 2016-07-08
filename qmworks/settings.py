
__all__ = ['Settings']

import plams

from noodles import Storable
from six import string_types


class Settings(plams.Settings, Storable):
    """
    This is a subclass of the :class:`plams.Settings`.
    TODO: explain the differences between this class and plams' Settings
    """

    def __getitem__(self, name):
        """Like regular __getitem__, but if name is a string, lowercase it."""
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            return dict.__getitem__(self, names[0]).__getitem__('.'.join(names[1:]))
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """Like regular __setitem__, but if name is a string, lowercase it."""
        if isinstance(value, dict):
            value = Settings(value)
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            dict.__getitem__(self, names[0]).__setitem__('.'.join(names[1:]), value)
        else:
            dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """Like regular __delitem__, but if name is a string, lowercase it."""
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            dict.__getitem__(self, names[0]).__delitem__('.'.join(names[1:]))
        else:
            dict.__delitem__(self, name)

    def copy(self):
        """
        """
        ret = Settings()
        for name in self:
            if isinstance(self[name], plams.Settings):
                ret[name] = self[name].copy()
            else:
                ret[name] = self[name]
        return ret

    def as_dict(self):
        return dict(self)

    @classmethod
    def from_dict(cls, **kwargs):
        return cls(**kwargs)

    def overlay(self, other):
        """
        Return new instance of |Settings| that is a copy of this instance
        updated with *other*.
        """
        ret = self.copy()
        ret.update(other)
        return ret

