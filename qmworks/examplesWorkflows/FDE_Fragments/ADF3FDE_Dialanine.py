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

settings = Settings ()
settings.functional = 'BP86'
settings.specific.adf.xc.accint = 2.0

settings.basis = 'SZ'
settings.specific.adf.basis.core = 'Large'

# supermolecule calculation
#supermol_results =  adf(templates.singlepoint.overlay(settings), rdkitTools.rdkit2plams(supermol), job_name = 'supermol_singlepoint')
supermol_results =  dftb(templates.singlepoint, rdkitTools.rdkit2plams(supermol), job_name = 'supermol_singlepoint')
print(run(supermol_results.dipole))
# supermol_dipole = supermol_results.get_dipole_vector()

#frags,caps = rdkitTools.partition_protein(supermol, cap=None)
frags, caps = rdopp.partition_protein(supermol, cap=None)

mfcc_result = run(mfcc(dftb, frags, caps, settings))

print(mfcc_result.get_dipole_vector())