# Default imports
from noodles import gather
from qmworks import (Settings, templates, run, molkit)

# User Defined imports
from qmworks import orca, dftb
from scm import plams
# ========== =============

template = "Cc1cc(C)cc(C)c1[PH+](c2c(C)cc(C)cc2C)C[BH3-]"

list_of_modifications = {"Me": "[#0:1]>>[CH3:1]"}
                         # "Ph": "[#0:1]>>[C:1]c1ccccc1",
                         # "tBu": "[#0:1]>>[C:1](C)(C)C",
                         # "PhF5": "[#0:1]>>[C:1]c1c(F)c(F)c(F)c(F)c1F",
                         # "Mes": "[#0:1]>>[C:1]c1c(C)cc(C)cc1C",
                         # "PhCF3_2": "[#0:1]>>[C:1]c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
                         # "Cl": "[#0:1]>>[Cl:1]",
                         # "CF3": "[#0:1]>>[C:1](F)(F)F"}

# HH = plams.Molecule("H_Mes2PCBH2_TS3series1.xyz")
HH = plams.Molecule("H_TS3.xyz")
HH.guess_bonds()
newmol = molkit.apply_template(HH, template)
# Change the 2 hydrogens to dummy atoms
temp = molkit.modify_atom(newmol, 47, '*')
temp = molkit.modify_atom(temp, 48, '*')
# Put coordinates of dummy atoms to 0.0, this will force generating new
# coordinates for the new atoms
temp.atoms[47].move_to((0.0, 0.0, 0.0))
temp.atoms[48].move_to((0.0, 0.0, 0.0))
# Replace dummy atoms by new groups, generate coordinates only for new atoms
job_list = []
for mod, smirks in list_of_modifications.items():
    print(mod)
    temp2, freeze = molkit.apply_reaction_smarts(temp, smirks, complete=True)
    # generate job
    s = Settings()
    s.freeze = freeze
    s.specific.dftb.dftb.resourcesdir = 'QUASINANO2015'
    partial_geometry = dftb(templates.geometry.overlay(s), temp2,
                            job_name="partial_opt_" + mod)
    freq_dftb = dftb(templates.freq.overlay(s), partial_geometry.molecule,
                     job_name="freq_" + mod)
    s = Settings()
    s.specific.orca.main = "B97-D"
    s.specific.orca.basis.basis = "_6_31G"
    s.specific.orca.basis.pol = "_d"
    s.inithess = freq_dftb.hessian
    ts = orca(templates.ts.overlay(s), freq_dftb.molecule,
              job_name="ts_" + mod)
    freq = orca(templates.freq.overlay(s), ts.molecule,
                job_name="freq" + mod)
    job_list.append(freq)

results = run(gather(*job_list), folder="Chris_results", n_processes=1)
