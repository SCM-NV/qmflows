# Default imports
from qmworks import (templates, run)
from qmworks import molkit
from qmworks.components import mfcc
from noodles import gather
from rdkit import Chem

# User Defined imports
from qmworks.packages.SCM import dftb

# from qmworks.components import adffragmentsjob

rdmol = Chem.MolFromPDBFile('Dialanine.pdb')
supermol = molkit.add_prot_Hs(rdmol)

# settings = Settings ()
# settings.functional = 'bp86'
#
# settings.basis = 'DZP'
# settings.specific.adf.basis.core = 'Large'

# supermolecule calculation
# supermol_job =  adf(templates.singlepoint.overlay(settings), supermol,
#                     job_name='supermol_singlepoint')
supermol_job = dftb(templates.singlepoint, supermol,
                    job_name='supermol_singlepoint')

frags, caps = molkit.partition_protein(supermol, cap=None)
# mfcc_job = mfcc(adf, frags, caps, settings)
mfcc_job = mfcc(dftb, frags, caps)

supermol_result, mfcc_result = run(gather(supermol_job, mfcc_job))

print(supermol_result.dipole)
print(mfcc_result.dipole)
