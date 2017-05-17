
__all__ = ['json2Settings']

# ================> Python Standard  and third-party <=========================
from qmworks.utils import dict2Setting

import json


def json2Settings(xs):
    """
    transform a string containing some data in JSON format
    to a Settings object
    """
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = json.loads(xs)  # Json object must be string
    return dict2Setting(s)
