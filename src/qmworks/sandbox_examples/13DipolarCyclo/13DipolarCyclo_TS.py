# Default imports
from qmworks import Settings, templates, run

# User Defined imports
from qmworks.packages.SCM import dftb
from qmworks.packages import orca
from qmworks.components import Distance
from scm import plams
import io

# ========== =============

# 5 membered ring from which ozone will dissociate
xyz_file = io.StringIO('''7

  C   -1.46639100895973     -0.66859400443971     -0.31255505682806
  O    0.00282072410259      0.01474573531737      1.30606888074830
  O    1.04388936077663     -0.01936067299255      0.53927051691306
  O    0.93874565776916      0.76907359603763     -0.48055185294869
  C   -0.93375384866586     -0.24065248892727     -1.32868508338709
  H   -2.06330249766107     -1.10524960979939      0.46177231972758
  H   -0.61425226052981      0.06238167064304     -2.30462642727095
''')

guessTS = plams.Molecule()
guessTS.readxyz(xyz_file, 1)

# Some dftb specific settings
dftb_set = Settings()
dftb_set.specific.dftb.dftb.scc

# Calculate the DFTB hessian
freq_dftb = dftb(templates.freq.overlay(dftb_set), guessTS, job_name="freq_dftb")

# User define Settings
settings = Settings()
settings.functional = "b3lyp"
settings.specific.orca.basis.basis = "_3_21G"
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
