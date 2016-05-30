
__author__ = "Felipe Zapata"

__all__ = ['configSSHKeys', 'check_rsa', 'generate_rsa_keys', 'copy_keys']

# ================> Python Standard  and third-party <==========

import fnmatch
import os
import os.path    as osp
import subprocess

# ===================> Environment <============================

home   = osp.expanduser("~")
sshdir = osp.join(home, ".ssh")

# =================> <====================


def configSSHKeys(hosts):
    """
    Check if there is a RSA key on the local host, otherwise create it and copy the
    credentials to the remote ``hosts``
    """
    check_rsa()
    for host in hosts:
        copy_keys(host)


def check_rsa():
    """
    Check if there is rsa key in the $HOME/.ssh folder,
    if there is no a key, generate one.
    """
    files = os.listdir(sshdir)
    rsa   = fnmatch.filter(files, "*rsa*")
    pub   = fnmatch.filter(files, "*rsa*pub")
    if rsa and pub:
        print("RSA key already exist on localhost")
        return pub[0]
    else:
        msg = 'The localhost does not have a RSA key in the folder:{}\nGenerating one'.format(sshdir)
        print(msg)
        return generate_rsa_keys()


def generate_rsa_keys():
    """
    Generates an 2048 RSA key and store it at $HOME/.ssh
    and return the name of the public key ``rsa.pub``
    """
    os.mkdir(sshdir, 384)  # 0600 in octal, in unix permission is set to 'drw-------'

    subprocess.call(["ssh-keygen", "-b", "2048", "-t", "rsa",
                     "-f", "$HOME/.ssh/rsa", "-q" "-N", ""])


def copy_keys(host):
    r = subprocess.call(["ssh-copy-id", host])
    if r != 0:
        err = 'can not copy RSA public key to host:{}\nManual intervention required'.format(host)
        raise RuntimeError(err)
