# Default imports
from qmworks import (Settings, templates, run, rdkitTools)
from noodles import gather

# User Defined imports
from math import sqrt
from qmworks.packages.SCM import (dftb, adf)
from qmworks.components import PES_scan, select_max

import plams
# ========== =============

plams.init()

hartree2kcal = 627.5095
startvalue = {'C': 2.1, 'N': 1.9, 'O': 1.8}

# User define Settings
settings = Settings()
settings.functional = "pbe"
settings.basis = "TZ2P"
settings.specific.dftb.dftb.scc


def calc_reactant(name):
    smiles = {'C': '[CH]#', 'N': '[N]#'}[name[0]] + "[N+]" +\
             {'C': '[CH2-]', 'N': '[NH-]', 'O': '[O-]'}[name[2]]
    mol = rdkitTools.smiles2plams(smiles)
    mol.properties.smiles = smiles
    reactant = adf(templates.geometry.overlay(settings),
                   dftb(templates.geometry.overlay(settings), mol,
                        job_name=name + "_reactant_DFTB").molecule,
                   job_name=name + "_reactant")
    return reactant


def calc_productAndTS(name):
    smiles = 'C1' + name[0] + "=" + name[1:] + 'C=1'
    mol = rdkitTools.smiles2plams(smiles)
    mol.properties.smiles = smiles
    product = adf(templates.geometry.overlay(settings),
                  dftb(templates.geometry.overlay(settings), mol,
                       job_name=name + "_product_DFTB").molecule,
                  job_name=name + "_product")

    constraint1 = "dist 1 2"
    constraint2 = "dist 4 5"

    sv1 = startvalue[name[0]]
    sv2 = startvalue[name[2]]

    # scan input
    if name[0] == name[2]:
        # symmetric molecule
        scan = {'constraint': [constraint1, constraint2],
                'surface': {'nsteps': 4, 'start': [sv1, sv2],
                            'stepsize': [0.1, 0.1]}}
    else:
        scan = {'constraint': constraint1,
                'surface': {'nsteps': 4, 'start': sv1, 'stepsize': 0.1},
                'scan': {'constraint': constraint2,
                         'surface': {'nsteps': 4, 'start': sv2, 'stepsize': 0.1}
                         }
                }

    # PES = gathered (promised) result objects for each point in the scan
    PES = PES_scan([dftb, adf], settings, product.molecule, scan,
                   job_name=name + "_PES")

    # get the object presenting the molecule with the maximum energy calculated from the scan
    apprTS = select_max(PES, 'energy')

    DFTB_freq = dftb(templates.freq.overlay(settings), apprTS.molecule,
                     job_name=name + "_DFTBfreq")

    t = Settings()
    t.specific.adf.geometry.inithess = DFTB_freq.archive.path

    # Run the TS optimization, using the default TS template
    TS = adf(templates.ts.overlay(settings).overlay(t), DFTB_freq.molecule,
             job_name=name + "_TS")

    adf_freq = adf(templates.freq.overlay(settings), TS.molecule,
                   job_name=name + "_freq")

    return product, adf_freq

# systematically generate list of names consisting of the 3 atoms that vary in
# the series of reactions
# e.g. "CNC", "CNN", etc.
reaction_names = set()
for a1 in ('C', 'N'):
    for a3 in ('C', 'N', 'O'):
        name = a1 + 'N' + a3
        reaction_names.add(name)
        print(name)

# generate all jobs
job_list = []
for name in reaction_names:
    reactant = calc_reactant(name)
    product, adf_freq = calc_productAndTS(name)
    job_list.append(gather(reactant, product, adf_freq))

# Need also the energy of ethene
acetylene = rdkitTools.smiles2plams('C#C')
E_acetylene_job = adf(templates.geometry.overlay(settings), acetylene,
                      job_name="acetylene").energy

# Actual execution of the jobs
E_acetylene, results = run(gather(E_acetylene_job, gather(*job_list)),
                           n_processes=1)


def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))

# extract table from results
table = {}
for reactant, product, TS in results:

    # Retrieve the molecular coordinates
    mol = TS.molecule
    d1 = bond_distance(mol.atoms[0].coords, mol.atoms[1].coords)
    d2 = bond_distance(mol.atoms[3].coords, mol.atoms[4].coords)

    Eact = (TS.energy - E_acetylene - reactant.energy) * hartree2kcal
    Ereact = (product.energy - E_acetylene - reactant.energy) * hartree2kcal
    smiles = reactant.molecule.properties.smiles
    table[smiles] = [Eact, Ereact, d1, d2]

# print table
print("Reactant                Eact  Ereact   Bond1   Bond2")
for smiles in table:
    print("{0:20s} {1:7.1f} {2:7.1f} {3:7.2f} \
    {4:7.2f}".format(smiles, *table[smiles]))
