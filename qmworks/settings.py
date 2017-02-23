
__all__ = ['Settings']

import plams

from noodles import Storable
from six import string_types


class Settings(plams.core.settings.Settings, Storable):
    """
    This is a subclass of the :class:`plams.core.settings.Settings`.
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
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            return dict.__getitem__(self, names[0]).__getitem__('.'.join(names[1:]))
        return dict.__getitem__(self, name)

    def __setitem__(self, name, value):
        """
        Like plams Settings.__setitem__, but
        "settings['a.b'] = 'c'" is equivalent to "settings['a']['b'] = 'c'"
        """
        if isinstance(value, dict):
            value = Settings(value)
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            dict.__getitem__(self, names[0]).__setitem__('.'.join(names[1:]), value)
        else:
            dict.__setitem__(self, name, value)

    def __delitem__(self, name):
        """
        Like plams Settings.__setitem__, but
        "del settings['a.b']" is equivalent to "del settings['a']['b'] = 'c'"
        """
        if isinstance(name, string_types) and name.find('.') > -1:
            names = name.split('.')
            dict.__getitem__(self, names[0]).__delitem__('.'.join(names[1:]))
        else:
            dict.__delitem__(self, name)

    def copy(self):
        ret = Settings()
        for name in self:
            if isinstance(self[name], plams.core.settings.Settings):
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

    def update(self, other):
        """
        Like PLAMS update, but:
        __block_replace = True results in removing all existing key value pairs
        in the current node of the settings tree

        For example:

        >>> from qmworks import Settings
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
        for name in other:
            if isinstance(other[name], Settings):
                if name not in self or not isinstance(self[name], Settings):
                    self[name] = other[name].copy()
                else:
                    # __block_replace can be used to remove all existing key value pairs in an input block
                    if ('__block_replace' in other[name] and other[name]['__block_replace'] == True) or \
                        ('__block_replace' in self[name] and self[name]['__block_replace'] == True):
                        self[name] = Settings()
                    self[name].update(other[name])
            else:
                self[name] = other[name]

