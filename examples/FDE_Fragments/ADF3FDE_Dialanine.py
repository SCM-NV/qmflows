# Default imports
from qmworks import (templates, run)
from qmworks import molkit
from qmworks.components import mfcc
from noodles import gather
from rdkit import Chem

# User Defined imports
from qmworks.packages.SCM import dftb

# For the purpose of the example, define the pdb block here.
# Could als be read from file.

dialanine_pdb = '''
HEADER    OXIDOREDUCTASE                          06-JUL-94   1PBE
TITLE     CRYSTAL STRUCTURE OF THE P-HYDROXYBENZOATE HYDROXYLASE-SUBSTRATE
TITLE    2 COMPLEX REFINED AT 1.9 ANGSTROMS RESOLUTION. ANALYSIS OF THE ENZYME-
TITLE    3 SUBSTRATE AND ENZYME-PRODUCT COMPLEXES
ATOM    937  N   ALA A 124      28.566 107.090  69.900  1.00 19.62           N
ATOM    938  CA  ALA A 124      28.940 106.368  71.140  1.00 19.55           C
ATOM    939  C   ALA A 124      30.169 105.591  70.700  1.00 21.08           C
ATOM    940  O   ALA A 124      29.935 104.701  69.852  1.00 22.09           O
ATOM    941  CB  ALA A 124      27.826 105.376  71.472  1.00 28.37           C
ATOM    942  N   ALA A 125      31.330 105.921  71.145  1.00 22.72           N
ATOM    943  CA  ALA A 125      32.573 105.296  70.663  1.00 22.82           C
ATOM    944  C   ALA A 125      33.073 104.318  71.705  1.00 26.31           C
ATOM    945  O   ALA A 125      32.743 104.415  72.876  1.00 22.69           O
ATOM    946  CB  ALA A 125      33.596 106.405  70.359  1.00 22.03           C
END
'''

rdmol = Chem.MolFromPDBBlock(dialanine_pdb)
# Could be read from file with
# rdmol = Chem.MolFromPDBFile('dialanine.pdb')

supermol = molkit.add_prot_Hs(rdmol)

# Calculate dipole normally
supermol_job = dftb(templates.singlepoint, supermol,
                    job_name='supermol_singlepoint')

# Calculate dipole with mfcc approach
frags, caps = molkit.partition_protein(supermol, cap=None)
mfcc_job = mfcc(dftb, frags, caps)

supermol_dipole, mfcc_dipole = run(gather(supermol_job.dipole, mfcc_job.dipole))

print(supermol_dipole)
print(mfcc_dipole)
