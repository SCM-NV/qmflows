
__all__ = ['freq', 'geometry', 'singlepoint', 'ts', 'get_template']

# ================> Python Standard  and third-party <==========

from os.path import join
import json
import pkg_resources as pkg
#  ==================> Internal Modules  <=====================
from qmworks.utils import dict2Setting
# ==================================================


def get_template(template_name):
    """
    The function :func:`pkg_resources.resource_string` search for the location
    of the ``template_name`` file inside the installation directory an returns
    a strings with its contain. This string is subsequently converted in a
    python *dict* object using :func:`json.loads` and finally recursively
    transform into a |Settings| object.
    Using this function there are created the default templates for the
    following common calculations:

    * Singlepoint
    * Geometry optimization
    * Transition State optimization
    * Frequencies
    """
    path = join("data/templates", template_name)
    xs = pkg.resource_string("qmworks", path)
    s = json.loads(xs.decode())  # Json object must be string
    return dict2Setting(s)


# Generic Templates
singlepoint = get_template('singlepoint.json')
geometry = get_template('geometry.json')
ts = get_template('ts.json')
freq = get_template('freq.json')
