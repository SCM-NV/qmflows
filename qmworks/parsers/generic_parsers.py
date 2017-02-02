
__all__ = ['awk_file', 'extract_line_value']

from pyparsing import  (OneOrMore, SkipTo, Suppress)
from qmworks.parsers.parser import parse_file

import os
import subprocess


def awk_file(filename, plams_dir=None, script='', progfile=None, **kwargs):
    """awk_file(filename, script='', progfile=None, **kwargs)
    Execute an AWK script on a file given by *filename*.

    The AWK script can be supplied in two ways: either by directly passing
    the contents of the script (should be a single string) as a *script*
    argument, or by providing the path (absolute or relative to the file
    pointed by *filename*) to some external file containing the actual AWK
    script using *progfile* argument. If *progfile* is not ``None``, the
    *script* argument is ignored.

    Other keyword arguments (*\*\*kwargs*) can be used to pass additional
    variables to AWK (see ``-v`` flag in AWK manual)

    Returned value is a list of lines (strings). See ``man awk`` for details.
    """
    cmd = ['awk']
    for k, v in kwargs.items():
        cmd += ['-v', '%s=%s' % (k, v)]
    if progfile:
        if os.path.isfile(progfile):
            cmd += ['-f', progfile]
        else:
            raise FileNotFoundError('File %s not present' % progfile)
    else:
        cmd += [script]

    new_cmd = cmd + [filename]
    ret = subprocess.check_output(new_cmd, cwd=plams_dir).decode('utf-8').split('\n')
    if ret[-1] == '':
        ret = ret[:-1]
    result = []
    for i in ret:
        try:
            v = int(i)
        except ValueError:
            try:
                v = float(i)
            except ValueError:
                v = i
        result.append(v)
    if len(result) == 1:
        result = result[0]
    return result


def extract_line_value(file_name, pattern=None, pos=0):
    """
    Get a field record from a file.
    Search for lines containing `pattern` and return the last line
    containing that value.

    :returns: value at position `pos` in the last found line, containing pattern.
    """
    parse_Line = OneOrMore(Suppress(SkipTo(pattern)) + SkipTo('\n'))
    properties = parse_file(parse_Line, file_name).asList()
    last_line = properties[-1].split()

    return float(last_line[pos])
