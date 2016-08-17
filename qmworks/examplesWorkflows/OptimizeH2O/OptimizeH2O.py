# Default imports
from qmworks import Settings, templates
from qmworks.draw_workflow import draw_workflow
from plams import Molecule
import plams

# User Defined imports
from qmworks.packages.SCM import adf, dftb
from qmworks.packages import run


plams.init()  # This is a plams requirement we like to get rid of
# global config
config.log.stdout = -1000

h2o = Molecule('h2o.xyz', 'xyz')

h2o_geometry = dftb(templates.geometry, h2o)

s = Settings()
# generic keyword "basis" must be present in the generic dictionary
s.basis = "DZP"
# "specific" allows the user to apply specific keywords for a
# package that are not in a generic dictionary
s.specific.adf.basis.core = "large"


h2o_singlepoint = adf(templates.singlepoint.overlay(s), h2o_geometry.molecule)

dipole = h2o_singlepoint.dipole

draw_workflow("wf.svg", dipole._workflow)

final_result = run(dipole, n_processes=1)

print(final_result)
