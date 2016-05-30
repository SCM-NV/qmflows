# import math
# import copy
# import numpy

from qmworks.quantumCodes.ADF import adf
from qmworks import Settings, pkg_templates
from plams import Molecule #, init, finish
from noodles import run_local
from qmworks.pkg_templates.inputGeneration import ADF_input_gen

### Plams needs those global config parameters
import builtins as btins
btins.config = Settings()
config.log.file = 0
config.log.stdout = 5
##############################################

supermol = Molecule('test_files/cys_cys.pdb', 'pdb')
supermol.spin = 0

us = Settings()
# us.functional = "BP86"   # these are generic keywords
us.basis = "SZ"
us.specific.adf.basis.core = "Large"    # these are ADF specific keyword
# us.specific.adf.xc.gga = "BP86"
us.specific.adf.integration.accint = 2.0

job_settings = pkg_templates.singlepoint.overlay(us)

# supermolecule calculation
supermol_results = adf(job_settings, supermol)
#supermol_dipole = supermol_results.dipole_vector

supermol_dipole = supermol_results.get_property('Dipole', 'Properties')

result = run_local.run(supermol_dipole)
print(result)

exit()

frags = partition_protein(supermol)

# MFCC calculation
mfcc_settings = Templates.mfcc.overlay(user_settings)
mfcc_results = Workflows.mfcc(frags, mfcc_settings, Adf)

mfcc_dipole = mfcc_results.dipole_vector
mfcc_frags  = mfcc_results.fragmentlist

# 3-FDE (0) calculation
# settings.set_ncycles(1000)
# settings.set_mixing(0.1)
# settings.set_convergence(converge=1.0e-4)
fde_settings = Templates.fde3.overlay(user_settings).overlay({
    "fde.RHO1FITTED": '', 
    "fde.CapDensConv": 1e-3,
    "fde.DENSTYPE": 'SCFfitted' })
fde_results = Workflows.fde(mfcc_frags, fde_settings, Adf)

results = run(fde_results, supermol_dipole)

fde_dipole = fde_results.dipole_vector
fde_dens= fde_results.get_density(spacing=0.2, fit=True)

supermol_dens = supermol_results.get_density(grid=fde_dens.grid, fit=True)

diff_dens = fde_dens - supermol_dens
diff_dens.get_cubfile(os.path.join(pyadfenv.outdir,'diff-dens.cub'))

print("Dipole moment (a.u.): ")
print("Supermolecule: ", supermol_dipole)
print("MFCC:          ", mfcc_dipole)
print("3-FDE:         ", fde_dipole)
print

sqrt_error = math.sqrt(diff_dens.integral(func=lambda x: x*x))
abs_error = diff_dens.integral(func=lambda x: abs(x))

print(" fde Diff Dens Int        : ", diff_dens.integral())
print(" fde Diff Dens Int (Sqr.) : ", sqrt_error)
print(" fde Diff Dens Int (Abs.) : ", abs_error)

testobj.assertAlmostEqual (sqrt_error, 0.00371608, 5)
testobj.assertAlmostEqual (abs_error,  0.04421127, 5)

def check_dipole (testobj, dipole_vect_ref, dipole_vect):

    testobj.assertAlmostEqual(dipole_vect[0], dipole_vect_ref[0], 4)
    testobj.assertAlmostEqual(dipole_vect[1], dipole_vect_ref[1], 4)
    testobj.assertAlmostEqual(dipole_vect[2], dipole_vect_ref[2], 4)

check_dipole (testobj, [-0.38075312, -0.45339314,  0.17361172], supermol_dipole)
check_dipole (testobj, [-0.37046423, -0.46913504,  0.19081185], mfcc_dipole)
check_dipole (testobj, [-0.39874326, -0.45753297,  0.211593841], fde_dipole)


