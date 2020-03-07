from qmflows import Settings, run, md, singlepoint, cp2k
from qmflows.packages.cp2k_mm import cp2k_mm

from scm.plams import readpdb


mol = readpdb('/Users/bvanbeek/Downloads/ff_assignment.028.pdb')
symbol_map = {at.properties.pdb_info.get('Name', at.symbol): at.symbol for at in mol}

s = Settings()
s.specific.cp2k += md.specific.cp2k_mm

job = cp2k_mm(mol=mol, settings=s)
results = run(job, path='/Users/bvanbeek/Downloads')
