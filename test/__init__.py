from .test_xyz_parser import *

import os

def exec_example(folder, filename):
    local_env = {}
    global_env = {}
    os.chdir(folder)
    try:
        exec(open(filename).read(), global_env, local_env)
        os.chdir('../..')
    except:
        # Make sure that if the tested code breaks the current dir is restored
        # Otherwise the subsequent tests will break for the wrong reason
        os.chdir('../..')
        raise RuntimeError()
    return local_env
