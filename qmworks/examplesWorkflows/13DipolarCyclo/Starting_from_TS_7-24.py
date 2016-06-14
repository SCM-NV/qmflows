# Default imports
from qmworks import Settings, templates, run, rdkitTools, draw_workflow
from noodles import gather
from plams import Molecule

# User Defined imports
from math import sqrt
from qmworks.packages.SCM import dftb,adf
#from qmworks.packages.SCM import dftb as adf  # This is for testing purposes
from qmworks.components import PES_scan, select_max

import plams
# ========== =============

plams.init()

cnc_ts = plams.Molecule('CNC_TS.mol')

# User define Settings
settings = Settings()
settings.functional = "pbe"
settings.basis = "TZ2P"
settings.specific.dftb.dftb.scc

job_list = []
for s1 in ('C','N','O'):
    n1=rdkitTools.modify_atom(cnc_ts,0,s1)
    for s2 in ('N','O','S'):
        n2=rdkitTools.modify_atom(n1,1,s2)
        for s3 in ('C','N','O'):
            n3=rdkitTools.modify_atom(n2,2,s3)
            name=s1+'-'+s2+'-'+s3
            n3.properties.name = name
            dftb_freq = dftb(templates.freq.overlay(settings),n3,job_name = name+"_dftb_freq")
            t = Settings()
            t.specific.adf.geometry.inithess = dftb_freq.archive.path
            job_list.append(adf(templates.ts.overlay(settings).overlay(t),n3,job_name = name+"_ts"))

wf = gather(*job_list)

# Actual execution of the jobs
results = run(wf, n_processes = 1)

def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))

# extract table from results
print("System                Bond1   Bond2")
for TS in results:

    # Retrieve the molecular coordinates
    mol = TS.molecule
    d1 = bond_distance(mol.atoms[0].coords, mol.atoms[4].coords)
    d2 = bond_distance(mol.atoms[3].coords, mol.atoms[2].coords)
    
    name = TS.molecule.properties.name
    print("{0:20s} {1:7.2f} {2:7.2f}".format(name,d1,d2))
