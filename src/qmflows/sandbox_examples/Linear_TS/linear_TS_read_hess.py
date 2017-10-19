# Default imports
from qmflows import Settings, templates, run
from noodles import gather
from plams import Molecule

# User Defined imports
from math import sqrt
from qmflows.packages.SCM import dftb, adf
from qmflows.components import PES_scan, select_max

import plams
# ========== =============


def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))
# ========== =============

plams.init()

# Read the Molecule from file
cnc = Molecule('C-N-C.mol', 'mol')

# User define Settings
settings = Settings()
settings.functional = "pbe"
settings.basis = "TZ2P"
settings.specific.dftb.dftb.scc

constraint1 = "dist 1 5"
constraint2 = "dist 3 4"

# scan input
scan = {'constraint': [constraint1, constraint2],
        'surface': {'nsteps': 6, 'start': [2.3, 2.3],
                    'stepsize': [0.1, 0.1]}}

# returns a set of results object containing the output of
# each point in the scan
lt = PES_scan([dftb,adf], settings, cnc, scan)

# Gets the object presenting the molecule
# with the maximum energy calculated from the scan
apprTS = select_max(lt, "energy")

appr_hess = dftb(templates.freq.overlay(settings), apprTS.molecule)

t = Settings()
t.specific.adf.geometry.inithess = appr_hess.archive.path

# Run the TS optimization with ADF, using initial hessian from DFTB freq calculation
ts = run(adf(templates.ts.overlay(settings).overlay(t), appr_hess.molecule), n_processes = 1)

# Retrieve the molecular coordinates
mol = ts.molecule
r1 = mol.atoms[0].coords
r2 = mol.atoms[4].coords

print("TS Bond distance:", bond_distance(r1, r2))
print("TS Energy:", ts.energy)

