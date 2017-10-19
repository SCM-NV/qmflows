# Default imports
from qmflows import (Settings, run)
from plams import Molecule

# User Defined imports
from math import sqrt
from qmflows.packages.SCM import dftb, adf
from qmflows.components import (Distance, PES, select_max)

# ========== =============


def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))
# ========== =============

# Read the Molecule from file
cnc = Molecule('C-N-C.mol', 'mol')


# User define Settings
settings = Settings()
settings.functional = "pbe"
settings.basis = "SZ"
settings.specific.dftb.dftb.scc

constraint1 = Distance(0, 4)
constraint2 = Distance(2, 3)

# scan input
pes = PES(cnc, constraints=[constraint1, constraint2],
          offset=[2.3, 2.3], get_current_values=False, nsteps=2, stepsize=[0.1, 0.1])

# returns a set of results object containing the output of
# each point in the scan
lt = pes.scan([dftb, adf], settings)

# Gets the object presenting the molecule
# with the maximum energy calculated from the scan
apprTS = select_max(lt, "energy")

# Run the TS optimization, using the default TS template
ts = run(apprTS)

# Retrieve the molecular coordinates
mol = ts.molecule
r1 = mol.atoms[0].coords
r2 = mol.atoms[4].coords

print("TS Bond distance:", bond_distance(r1, r2))
print("TS Energy:", ts.energy)
