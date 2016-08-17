
__all__ = ['json2Settings', 'search_environ_var']

# ================> Python Standard  and third-party <=========================
from qmworks.utils import dict2Setting

import json
import os
# ====================================<>=======================================


def json2Settings(xs):
    """
    transform a string containing some data in JSON format
    to a Settings object
    """
    if isinstance(xs, bytes):
        xs = xs.decode()
    s = json.loads(xs)  # Json object must be string
    return dict2Setting(s)

# ==========================> Files Utilities <================================


def search_environ_var(var):
    """
    Looks if the environmental variable ``var`` is defined.
    """
    try:
        defaultPath = os.environ[var]
    except KeyError:
        pass

    if not defaultPath:
        msg = 'There is not an environmental variable called: {}'.format(var)
        raise EnvironmentError(msg)

    return defaultPath

