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

supermol_job = dftb(templates.singlepoint, supermol,
                    job_name='supermol_singlepoint')

frags, caps = molkit.partition_protein(supermol, cap=None)
mfcc_job = mfcc(dftb, frags, caps)

supermol_dipole, mfcc_dipole = run(gather(supermol_job.dipole, mfcc_job.dipole))

print(supermol_dipole)
print(mfcc_dipole)
