# Default imports
from qmworks import Settings, templates, run, rdkitTools
from qmworks import rdkitTools as rdopp
from qmworks.components import mfcc
from noodles import gather, schedule, Storable
import plams
from rdkit import Chem

plams.init()

# User Defined imports
from qmworks.packages.SCM import dftb, adf

# from qmworks.components import adffragmentsjob

rdmol = Chem.MolFromPDBFile('Dialanine.pdb')
supermol = rdopp.add_prot_Hs(rdmol)

# settings = Settings ()
# settings.functional = 'bp86'
#
# settings.basis = 'DZP'
# settings.specific.adf.basis.core = 'Large'

# supermolecule calculation
# supermol_job =  adf(templates.singlepoint.overlay(settings), supermol, job_name = 'supermol_singlepoint')
supermol_job =  dftb(templates.singlepoint, supermol, job_name = 'supermol_singlepoint')
# supermol_dipole = supermol_results.get_dipole_vector()

#frags,caps = rdkitTools.partition_protein(supermol, cap=None)
frags, caps = rdopp.partition_protein(supermol, cap=None)
# mfcc_job = mfcc(adf, frags, caps, settings)
mfcc_job = mfcc(dftb, frags, caps)

supermol_result, mfcc_result = run(gather(supermol_job, mfcc_job))

print(supermol_result.dipole)
print(mfcc_result.dipole)