# Default imports
from qmworks import Settings, templates, run, rdkitTools
from qmworks.draw_workflow import draw_workflow
from noodles import gather
# from plams import Molecule

# User Defined imports
from math import sqrt
from qmworks.packages.SCM import dftb, adf
# from qmworks.packages.SCM import dftb as adf  # This is for testing purposes
from qmworks.components import Distance, PES, select_max

import plams
# ========== =============

plams.init()
config.log.stdout = -1000  # noqa

hartree2kcal = 627.5095

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
    reactant = adf(
        templates.geometry.overlay(settings),
        dftb(templates.geometry.overlay(settings), mol,
             job_name=name + "_reactant_DFTB").molecule,
        job_name=name + "_reactant")
    return reactant


def calc_productAndTS(name):
    smiles = 'C1' + name[0] + "=" + name[1:] + 'C1'
    mol = rdkitTools.smiles2plams(smiles)
    mol.properties.smiles = smiles
    product = adf(
        templates.geometry.overlay(settings),
        dftb(templates.geometry.overlay(settings),
             mol, job_name=name + "_product_DFTB").molecule,
        job_name=name + "_product")

    constraint1 = Distance(1, 0)
    constraint2 = Distance(4, 3)

    # scan input
    pes = PES(product.molecule, constraints=[constraint1, constraint2],
              offset=[2.0, 2.0], get_current_values=False, nsteps=5, stepsize=[0.1, 0.1])

    # PES = gathered (promised) result objects for each point in the scan
    pesscan = pes.scan([dftb, adf], settings, job_name=name + "_PES")

    # get the object presenting the molecule with the maximum energy
    # calculated from the scan
    apprTS = select_max(pesscan, 'energy')

    DFTB_freq = dftb(
        templates.freq.overlay(settings), apprTS.molecule,
        job_name=name + "_DFTBfreq")

    t = Settings()
    t.specific.adf.geometry.inithess = DFTB_freq.kf.path

    # Run the TS optimization, using the default TS template
    TS = adf(
        templates.ts.overlay(settings).overlay(t),
        DFTB_freq.molecule, job_name=name + "_TS")

    adf_freq = adf(
        templates.freq.overlay(settings), TS.molecule,
        job_name=name + "_freq")

    return product, adf_freq

# systematically generate list of names consisting of the 3 atoms
# that vary in the series of reactions
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
ethene = rdkitTools.smiles2rdkit('C=C')
E_ethene_job = adf(
    templates.geometry.overlay(settings),
    ethene, job_name="ethene").energy

# finalize and draw workflow
wf = gather(E_ethene_job, gather(*job_list))
draw_workflow("wf.svg", wf._workflow)

# Actual execution of the jobs
E_ethene, results = run(wf, n_processes=2)


def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))

# extract table from results
table = {}
for reactant, product, TS in results:
    # Retrieve the molecular coordinates
    mol = TS.molecule
    d1 = bond_distance(mol.atoms[0].coords, mol.atoms[1].coords)
    d2 = bond_distance(mol.atoms[3].coords, mol.atoms[4].coords)

    Eact = (TS.energy - E_ethene - reactant.energy) * hartree2kcal
    Ereact = (product.energy - E_ethene - reactant.energy) * hartree2kcal
    smiles = reactant.molecule.properties.smiles
    table[smiles] = [Eact, Ereact, d1, d2]

# print table
print("Reactant                Eact  Ereact   Bond1   Bond2")
for smiles in table:
    print("{0:20s} {1:7.1f} {2:7.1f} {3:7.2f} {4:7.2f}".format(
        smiles, *table[smiles]))
