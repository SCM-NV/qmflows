# Default imports
from qmworks import Settings, templates, run, molkit
from qmworks.draw_workflow import draw_workflow
from noodles import gather

# User Defined imports
from qmworks.packages.SCM import dftb
from qmworks.packages.SCM import dftb as adf  # This is for testing purposes
from qmworks.components import Distance, PES, select_max

from scm import plams
# ========== =============

plams.init()
config.log.stdout = -1000  # noqa

hartree2kcal = 627.5095

# List of reactions, defined by name, reactants and products (as smiles strings)
reactions = [
    ["et01", "C=C", "[CH]#[N+][CH2-]", "C1C=NCC1"],
    ["et02", "C=C", "[CH]#[N+][NH-]", "C1C=NNC1"],
    ["et03", "C=C", "[CH]#[N+][O-]", "C1C=NOC1"],
    ["et04", "C=C", "[N]#[N+][CH2-]", "C1N=NCC1"],
    ["et05", "C=C", "[N]#[N+][NH-]", "C1N=NNC1"],
    ["et06", "C=C", "[N]#[N+][O-]", "C1N=NOC1"],
    ["et07", "C=C", "[CH2]=[NH+][CH2-]", "C1CNCC1"],
    ["et08", "C=C", "[CH2]=[NH+][NH-]", "C1CNNC1"],
    ["et09", "C=C", "[CH2]=[NH+][O-]", "C1CNOC1"],
    ["et10", "C=C", "[NH]=[NH+][NH-]", "C1NNNC1"],
    ["et11", "C=C", "[NH]=[NH+][O-]", "C1NNOC1"],
    ["et12", "C=C", "O=[NH+][O-]", "C1ONOC1"],
    ["et13", "C=C", "[CH2]=[O+][CH2-]", "C1COCC1"],
    ["et14", "C=C", "[CH2]=[O+][NH-]", "C1CONC1"],
    ["et15", "C=C", "[CH2]=[O+][O-]", "C1COOC1"],
    ["et16", "C=C", "[NH]=[O+][NH-]", "C1NONC1"],
    ["et17", "C=C", "[NH]=[O+][O-]", "C1NOOC1"],
    ["et18", "C=C", "O=[O+][O-]", "C1OOOC1"],
    ["et19", "C=C", "[CH2]=[S+][CH2-]", "C1CSCC1"],
    ["et20", "C=C", "[CH2]=[S+][NH-]", "C1CSNC1"],
    ["et21", "C=C", "[CH2]=[S+][O-]", "C1CSOC1"],
    ["et22", "C=C", "[NH]=[S+][NH-]", "C1NSNC1"],
    ["et23", "C=C", "[NH]=[S+][O-]", "C1NSOC1"],
    ["et24", "C=C", "O=[S+][O-]", "C1OSOC1"],

    ["ac01", "C#C", "[CH]#[N+][CH2-]", "C1C=NCC=1"],
    ["ac02", "C#C", "[CH]#[N+][NH-]", "C1C=NNC=1"],
    ["ac03", "C#C", "[CH]#[N+][O-]", "C1C=NOC=1"],
    ["ac04", "C#C", "[N]#[N+][CH2-]", "C1N=NCC=1"],
    ["ac05", "C#C", "[N]#[N+][NH-]", "C1N=NNC=1"],
    ["ac06", "C#C", "[N]#[N+][O-]", "C1N=NOC=1"],
    ["ac07", "C#C", "[CH2]=[NH+][CH2-]", "C1CNCC=1"],
    ["ac08", "C#C", "[CH2]=[NH+][NH-]", "C1CNNC=1"],
    ["ac09", "C#C", "[CH2]=[NH+][O-]", "C1CNOC=1"],
    ["ac10", "C#C", "[NH]=[NH+][NH-]", "C1NNNC=1"],
    ["ac11", "C#C", "[NH]=[NH+][O-]", "C1NNOC=1"],
    ["ac12", "C#C", "O=[NH+][O-]", "C1ONOC=1"],
    ["ac13", "C#C", "[CH2]=[O+][CH2-]", "C1COCC=1"],
    ["ac14", "C#C", "[CH2]=[O+][NH-]", "C1CONC=1"],
    ["ac15", "C#C", "[CH2]=[O+][O-]", "C1COOC=1"],
    ["ac16", "C#C", "[NH]=[O+][NH-]", "C1NONC=1"],
    ["ac17", "C#C", "[NH]=[O+][O-]", "C1NOOC=1"],
    ["ac18", "C#C", "O=[O+][O-]", "C1OOOC=1"],
    ["ac19", "C#C", "[CH2]=[S+][CH2-]", "C1CSCC=1"],
    ["ac20", "C#C", "[CH2]=[S+][NH-]", "C1CSNC=1"],
    ["ac21", "C#C", "[CH2]=[S+][O-]", "C1CSOC=1"],
    ["ac22", "C#C", "[NH]=[S+][NH-]", "C1NSNC=1"],
    ["ac23", "C#C", "[NH]=[S+][O-]", "C1NSOC=1"],
    ["ac24", "C#C", "O=[S+][O-]", "C1OSOC=1"]]

# Bonds to be formed (based on atoms numbers in de product)
bond1 = Distance(1, 0)
bond2 = Distance(4, 3)

# User define Settings
settings = Settings()
settings.functional = "pbe"
settings.basis = "TZ2P"
settings.specific.dftb.dftb.scc


job_list = []
# Loop over all reactions
for name, reactant1, reactant2, product in reactions[:1]:

  # Prepare reactant1 job
    r1mol = molkit.from_smiles(reactant1)
    r1job = adf(
        templates.geometry.overlay(settings),
        dftb(templates.geometry.overlay(settings), r1mol,
             job_name=name + "_r1_DFTB").molecule,
        job_name=name + "_r1")

  # Prepare reactant2 job
    r2mol = molkit.from_smiles(reactant2)
    r2job = adf(
        templates.geometry.overlay(settings),
        dftb(templates.geometry.overlay(settings), r2mol,
             job_name=name + "_r2_DFTB").molecule,
        job_name=name + "_r2")

  # Prepare product job
    pmol = molkit.from_smiles(product)
    pmol.properties.name = name
    pjob = adf(
        templates.geometry.overlay(settings),
        dftb(templates.geometry.overlay(settings),
             pmol, job_name=name + "_p_DFTB").molecule,
        job_name=name + "_p")

  # Prepare TS jobs
    # Define PES
    pes = PES(pjob.molecule, constraints=[bond1, bond2],
              offset=[2.0, 2.0], get_current_values=False, nsteps=5, stepsize=[0.1, 0.1])

    # pesscan: gathered (promised) result objects for each point in the scan
    pesscan = pes.scan([dftb, adf], settings, job_name=name + "_PES")

    # Get the result with the maximum energy
    apprTS = select_max(pesscan, 'energy')

    # Calculate the DFTB hessian
    DFTBfreq = dftb(templates.freq.overlay(settings), apprTS.molecule, job_name=name + "_DFTBfreq")

    # Run the TS optimization, using the initial hession from DFTB
    t = Settings()
    t.inithess = DFTBfreq.hessian
    TS = adf(templates.ts.overlay(settings).overlay(t), DFTBfreq.molecule, job_name=name + "_TS")

    # Perform a freq calculation
    TSfreq = adf(templates.freq.overlay(settings), TS.molecule, job_name=name + "_freq")

  # Add the jobs to the job list
    job_list.append(gather(r1job, r2job, pjob, TSfreq, TS.optcycles))

# Finalize and draw workflow
wf = gather(*job_list)
draw_workflow("wf.svg", wf._workflow)

# Actual execution of the jobs
results = run(wf, n_processes=1)

# Extract table from results
table = {}
for r1result, r2result, presult, tsresult, optcycles in results:
    # Retrieve the molecular coordinates
    mol = tsresult.molecule
    d1 = bond1.get_current_value(mol)
    d2 = bond2.get_current_value(mol)

    Eact = (tsresult.energy - r1result.energy - r2result.energy) * hartree2kcal
    Ereact = (presult.energy - r1result.energy - r2result.energy) * hartree2kcal
    name = presult.molecule.properties.name
    smiles = presult.molecule.properties.smiles
    nnegfreq = sum([f < 0 for f in tsresult.frequencies])
    table[name] = [smiles, Eact, Ereact, d1, d2, nnegfreq, optcycles]

# Print table
print("Reaction Productsmiles    Eact  Ereact   Bond1   Bond2 NNegFreq OptCycles")
for name in sorted(table):
    print("{0:8s} {1:13s} {2:7.1f} {3:7.1f} {4:7.2f} {5:7.2f} {6:8d} {7:8d}".format(
        name, *table[name]))
