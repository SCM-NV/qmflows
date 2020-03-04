
__all__ = ['json2Settings', 'yaml2Settings']

import json
import yaml

from qmflows.utils import dict2Setting


def json2Settings(xs):
    """Transform a string containing some data in JSON format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = json.loads(xs)  # Json object must be string
    return dict2Setting(s)


def yaml2Settings(xs):
    """Transform a string containing some data in .yaml format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = yaml.load(xs, Loader=yaml.FullLoader)  # yaml object must be string
    return dict2Setting(s)
