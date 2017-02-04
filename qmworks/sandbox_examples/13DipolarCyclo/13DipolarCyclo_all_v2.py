# Default imports
from qmworks import Settings, templates, run, molkit
from qmworks.draw_workflow import draw_workflow
from noodles import gather

# User Defined imports
from qmworks.packages.SCM import dftb, adf
from qmworks.packages.orca import orca as adf
# from qmworks.packages.SCM import dftb as adf  # This is for testing purposes
from qmworks.components import Distance, select_max, select_min

import plams
# ========== =============

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

# Bonds to be formed (based on atom numbers in de product)
bond1 = Distance(0, 1)
bond2 = Distance(3, 4)

# User define Settings
settings = Settings()
settings.functional = "b3lyp"
settings.specific.orca.basis.basis = "_6_31G"
settings.specific.orca.basis.pol = "_d"
settings.specific.orca.basis.diff = "_p"
settings.specific.dftb.dftb.scc.ndiis = 4
settings.specific.dftb.dftb.scc.Mixing = 0.1
settings.specific.dftb.dftb.scc.iterations = 300


job_list = []
# Loop over all reactions
for name, r1_smiles, r2_smiles, p_smiles in reactions:

  # Prepare reactant1 job
    r1_mol =  molkit.from_smiles(r1_smiles)
    r1_dftb = dftb(templates.geometry.overlay(settings), r1_mol, job_name=name + "_r1_DFTB")
    r1 =      adf(templates.geometry.overlay(settings), r1_dftb.molecule, job_name=name + "_r1")
    r1_freq = adf(templates.freq.overlay(settings), r1.molecule, job_name=name + "_r1_freq")

  # Prepare reactant2 job
    r2s_mol =  molkit.from_smiles(r2_smiles, nconfs=10, forcefield='mmff', rms=0.1)
    r2s_dftb = [dftb(templates.geometry.overlay(settings), r2_mol, job_name=name + "_r2_DFTB") for r2_mol in r2s_mol]
    r2_dftb = select_min(gather(*r2s_dftb), 'energy')
    r2 =      adf(templates.geometry.overlay(settings), r2_dftb.molecule, job_name=name + "_r2")
    r2_freq = adf(templates.freq.overlay(settings), r2.molecule, job_name=name + "_r2_freq")

# Prepare product job
    ps_mol =  molkit.from_smiles(p_smiles, nconfs=10, forcefield='mmff', name=name, rms=0.1)
    ps_dftb = [dftb(templates.geometry.overlay(settings), p_mol, job_name=name + "_ps_DFTB") for p_mol in ps_mol]
    p_dftb = select_min(gather(*ps_dftb), 'energy')
    p =      adf(templates.geometry.overlay(settings), p_dftb.molecule, job_name=name + "_p")
    p_freq = adf(templates.freq.overlay(settings), p.molecule, job_name=name + "_p_freq")

# Prepare scan
    pes_jobs = []
    for d in range(6):
        consset = Settings()
        consset.constraint.update(bond1.get_settings(2.0 + d * 0.1))
        consset.constraint.update(bond2.get_settings(2.0 + d * 0.1))

        pes_name = name + "_pes_" + str(d)
        pes_dftb = dftb(templates.geometry.overlay(settings).overlay(consset), p.molecule, job_name=pes_name + "_DFTB")
        pes =      adf(templates.singlepoint.overlay(settings), pes_dftb.molecule, job_name=pes_name)
        pes_jobs.append(pes)

  # Get the result with the maximum energy
    apprTS = select_max(gather(*pes_jobs), 'energy')
  # Calculate the DFTB hessian
    DFTBfreq = dftb(templates.freq.overlay(settings), apprTS.molecule, job_name=name + "_freq_DFTB")
  # Run the TS optimization, using the initial hession from DFTB
    t = Settings()
    t.inithess = DFTBfreq.hessian
    TS = adf(templates.ts.overlay(settings).overlay(t), DFTBfreq.molecule, job_name=name + "_TS")

  # Perform a freq calculation
  #   freq_setting = Settings()
  #   freq_setting.specific.adf.geometry.frequencies = ""
  #   freq_setting.specific.adf.geometry.__block_replace = True
  #   TSfreq = adf(templates.geometry.overlay(settings).overlay(freq_setting), TS.molecule, job_name=name + "_freq")
    TSfreq = adf(templates.freq.overlay(settings), TS.molecule, job_name=name + "_freq")

  # Add the jobs to the job list
    job_list.append(gather(r1_freq, r2_freq, p_freq, TS, TSfreq))

# Finalize and draw workflow
wf = gather(*job_list)
draw_workflow("wf.svg", wf._workflow)

# Actual execution of the jobs
results = run(wf, n_processes=1)

# Extract table from results
table = {}
for r1_result, r2_result, p_result, ts_opt, ts_result in results:
    # Retrieve the molecular coordinates
    mol = ts_opt.molecule
    d1 = bond1.get_current_value(mol)
    d2 = bond2.get_current_value(mol)

    Eact = (ts_result.enthalpy - r1_result.enthalpy - r2_result.enthalpy) * hartree2kcal
    Ereact = (p_result.enthalpy - r1_result.enthalpy - r2_result.enthalpy) * hartree2kcal
    name = p_result.molecule.properties.name
    smiles = p_result.molecule.properties.smiles
    nimfreq = sum([f < 0 for f in ts_result.frequencies])
    table[name] = [smiles, Eact, Ereact, d1, d2, nimfreq, ts_opt.optcycles, ts_opt.runtime]

# Print table
print("Reaction Productsmiles    Eact  Ereact   Bond1   Bond2 NNegFreq TSoptCycles TSoptTime")
for name in sorted(table):
    print("{0:8s} {1:13s} {2:7.1f} {3:7.1f} {4:7.2f} {5:7.2f} {6:8d} {7:10d} {8:10.2f}".format(
        name, *table[name]))
