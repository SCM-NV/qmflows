__author__ = "Felipe Zapata"

__all__ = ['portForward']

# ================> Python Standard  and third-party <==========
from sys import platform as _platform
import subprocess


# =================> Internal Modules <================================

from qmworks.utils import flatten

# ======================================================================


def portForward(hostName, port='27017'):
    """
    Forward the traffic from ``port`` to ``port`` at remote ``hostname``.
    In unix the command: ssh -L 27017:localhost:27017 hostname -N
    is executed in a new process.

    """
    hostport   = '{}:localhost:{}'.format(port, port)
    
    unixSSH    = ['ssh', '-L', hostport]
    windowsSSH = []
    sshCommand = unixSSH if any([_platform == x
                                 for x in ["linux", "linux2", "darwin"]]) else windowsSSH
    cmd        = flatten([sshCommand, [hostName], ['-N']])

    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

