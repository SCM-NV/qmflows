import shutil
import pathlib
from os.path import isdir

from qmflows.utils import init_restart, InitRestart
from scm.plams import init, finish

PATH = pathlib.Path('test') / 'test_files'


def test_restart_init() -> None:
    """Tests for :func:`restart_init` and :class:`RestartInit`."""
    workdir = PATH / 'plams_workdir'
    try:
        init(PATH)
        finish()
        assert isdir(workdir)

        init_restart(PATH)
        assert isdir(workdir)
        assert not isdir(f'{workdir}.002')

        with InitRestart(PATH):
            assert isdir(workdir)
            assert not isdir(f'{workdir}.002')

    finally:
        shutil.rmtree(workdir) if isdir(workdir) else None
        shutil.rmtree(f'{workdir}.002') if isdir(f'{workdir}.002') else None
