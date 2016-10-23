# Default imports
from qmworks import (templates, run)
from qmworks import molkit
from qmworks.components import mfcc
from noodles import gather
from rdkit import Chem

# User Defined imports
from qmworks.packages.SCM import dftb

# from qmworks.components import adffragmentsjob

# It turned out important to rename the terminal carboxyl oxygens
# to which hydrogens were connected to "OXT" in order for RDKIT
# to correctly interpret the connectivity of the cys_cys.pdb
rdmol = Chem.MolFromPDBFile('cys_cys.pdb', removeHs=False)
supermol = rdmol

# Calculate dipole normally
supermol_job = dftb(templates.singlepoint, supermol,
                    job_name='supermol_singlepoint')
supermol_dipole = supermol_job.dipole

frags, caps = molkit.partition_protein(supermol, cap=None)
mfcc_job = mfcc(dftb, frags, caps)

supermol_dipole, mfcc_dipole = run(gather(supermol_job.dipole, mfcc_job.dipole))

print(supermol_dipole)
print(mfcc_dipole)
