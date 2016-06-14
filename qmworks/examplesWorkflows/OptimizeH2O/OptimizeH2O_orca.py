# Default imports
from qmworks import Settings, templates
from plams import Molecule
import plams

# User Defined imports
from qmworks.packages.SCM import adf, dftb
from qmworks.packages.orca import orca
from qmworks.packages import run


plams.init()  # This is a plams requirement we like to get rid of

h2o = Molecule('h2o.xyz', 'xyz', charge=0, multiplicity=1)

h2o_geometry = dftb(templates.geometry, h2o)

s = Settings()
# generic keyword "basis" must be present in the generic dictionary
s.basis = "sto_dzp"
# "specific" allows the user to apply specific keywords for a
# package that are not in a generic dictionary
# s.specific.adf.basis.core = "large"


h2o_singlepoint = orca(templates.singlepoint.overlay(s), h2o_geometry.molecule)

dipole = h2o_singlepoint.dipole

final_result = run(dipole, n_processes=1)

print(final_result)
