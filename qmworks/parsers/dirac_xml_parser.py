

import xml.etree.ElementTree as ET


def get_text(xs, read=float):
    """
    Read a Numerical value from the Leaf of a tree
    """
    if isinstance(xs, list):
        rs = list(map(lambda x: read(getattr(x, 'text')), xs))
    else:
        rs = read(xs.text)

    return rs


def parse_xml(filename, section=None, select_value=None):
    """
    Get a Property from the dirac xml output.
    :param filename: Name of the output file
    :param prop: Property to look up in the ouput

    :returns: Any
    """
    # whole output
    tree = ET.parse(filename)
    root = tree.getroot()

    # Property value
    xs = root.findall(section)

    # If there is a list of values and only the last is important
    if isinstance(xs, list) and select_value == 'last':
        xs = xs[-1]

    return get_text(xs)
