"""Use the AWK command line to read output files."""

__all__ = ['awk_file', 'extract_line_value', 'extract_line_values']

import subprocess
from typing import Union, List, Optional, Any

import numpy as np
from pyparsing import OneOrMore, SkipTo, Suppress

from .parser import parse_file
from ..type_hints import PathLike


Scalar = Union[int, float, str]


def awk_file(filename: PathLike, script: str = '',
             progfile: Optional[PathLike] = None
             ) -> Union[Scalar, List[Scalar]]:
    r"""Execute an AWK script on a file given by *filename*.

    The AWK script can be supplied in two ways: either by directly passing
    the contents of the script (should be a single string) as a *script*
    argument, or by providing the path (absolute or relative to the file
    pointed by *filename*) to some external file containing the actual AWK
    script using *progfile* argument. If *progfile* is not ``None``, the
    *script* argument is ignored.

    Returned value is a list of lines (strings). See ``man awk`` for details.
    """
    if not isinstance(filename, str):
        filename = str(filename)

    cmd = ['awk', script, filename]
    # if progfile:
    #     if os.path.isfile(progfile):
    #         cmd += ['-f', str(progfile)]
    #     else:
    #         raise FileNotFoundError('File %s not present' % progfile)
    # else:
    #     cmd += [script]

    # new_cmd = cmd + [filename]
    ret = subprocess.check_output(cmd).decode('utf-8').split('\n')
    if ret[-1] == '':
        ret = ret[:-1]

    result: List[Scalar] = []
    for i in ret:
        try:
            v: Scalar = int(i)
        except ValueError:
            try:
                v = float(i)
            except ValueError:
                v = i
        result.append(v)
    return result[0] if len(result) == 1 else result


def extract_line_value(file_name: PathLike,
                       pattern: Optional[str] = None,
                       pos: int = 0) -> float:
    """Get a field record from a file.

    Search for lines containing `pattern` and return the last line
    containing that value.
    :returns: value at position `pos` in the last found line, containing pattern.
    """
    parse_Line = OneOrMore(Suppress(SkipTo(pattern)) + SkipTo('\n'))
    properties = parse_file(parse_Line, file_name).asList()
    last_line = properties[-1].split()

    return float(last_line[pos])


def extract_line_values(
    file_name: PathLike,
    pattern: Optional[str] = None,
    pos: int = 0,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    dtype: Any = np.float64,
) -> np.ndarray:
    """Get multiply field records from a file.

    Search for lines containing `pattern` and return all line
    containing that value in the range defined by `start`, `stop` and `step`.
    :returns: value at position `pos` in all found lines, containing pattern.
    """
    parse_Line = OneOrMore(Suppress(SkipTo(pattern)) + SkipTo('\n'))
    properties = parse_file(parse_Line, file_name)
    iterator = (i.split()[pos] for i in properties[start:stop:step])
    return np.fromiter(iterator, dtype=dtype)
