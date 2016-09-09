

from pyparsing import (Literal, Regex, SkipTo, Suppress)


def parser_gamess(path_to_dat, prop):
    """
    Read numerical data from the `*.dat` file.
    """
    parser = {'totalenergy': parse_energy}

    try:
        return parse_file(parser[prop], path_to_dat)
    except KeyError:
        msg = "There is not parser from propiety: {}".format(prop)
        print(msg)
        raise


def parse_file(parser_fun, path_to_dat):
    return parser_fun().parseFile(path_to_dat)


def parse_energy():
    """
    parse Total energy.
    """
    Suppress(SkipTo(Regex(r'^E'))
