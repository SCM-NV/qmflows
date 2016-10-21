# Default imports
from qmworks import Settings, templates, rdkitTools
from qmworks.draw_workflow import draw_workflow
from plams import Molecule
import plams

from noodles import gather

# User Defined imports
from qmworks.packages.SCM import adf, dftb
from qmworks.packages.orca import orca
from qmworks.packages.gamess import gamess
from qmworks.packages import run

import numpy as np

plams.init()  # This is a plams requirement we like to get rid of

h2o = rdkitTools.smiles2plams('O')
h2o.properties.symmetry = 'C1'

h2o_freq = dftb(templates.freq, h2o, job_name = "freq").hessian
# h2o_freq = gamess(templates.freq, h2o, job_name = "freq", work_dir = "/home/lars/scr").hessian

s = Settings()
s.inithess = h2o_freq

h2o_opt = orca(templates.geometry.overlay(s), h2o, job_name = "opt")


energy = h2o_opt.energy
dipole = h2o_opt.dipole


wf = gather(energy, dipole)
draw_workflow('wf.svg', wf._workflow)

print(run(wf))


