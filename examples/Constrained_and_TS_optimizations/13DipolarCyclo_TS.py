# Default imports
from qmworks import Settings, templates, run

# User Defined imports
from qmworks.packages.SCM import dftb
from qmworks.packages import orca
from qmworks.components import Distance

import plams
# ========== =============

# 5 membered ring from which ozone will dissociate
guessTS = plams.Molecule('ac18_ts.xyz')

# Some dftb specific settings
dftb_set = Settings()
dftb_set.specific.dftb.dftb.scc

# Calculate the DFTB hessian
freq_dftb = dftb(templates.freq.overlay(dftb_set), guessTS, job_name="freq_dftb")

# User define Settings
settings = Settings()
settings.functional = "b3lyp"
settings.specific.orca.basis.basis = "_6_31G"
settings.specific.orca.basis.pol = "_d"
settings.specific.orca.basis.diff = "_p"
# Run the TS optimization, using the initial hessian from DFTB
settings.inithess = freq_dftb.hessian

# Define job
ts = orca(templates.ts.overlay(settings), guessTS, job_name="ts")

# Actual execution of the jobs
ts_result = run(ts)

# Bonds to be formed (based on atom numbers in de product)
bond1 = Distance(0, 1)
bond2 = Distance(3, 4)
d1 = bond1.get_current_value(ts_result.molecule)
d2 = bond2.get_current_value(ts_result.molecule)

optcycles = ts_result.optcycles

print('Distance 1 in ts:', d1)
print('Distance 2 in ts:', d2)
print('Optimization cycles:', optcycles)
