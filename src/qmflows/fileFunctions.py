
__all__ = ['json2Settings']

import json

from qmflows.utils import dict2Setting


def json2Settings(xs):
    """Transform a string containing some data in JSON format to a Settings object."""
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = json.loads(xs)  # Json object must be string
    return dict2Setting(s)
