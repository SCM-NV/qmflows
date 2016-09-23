# Default imports
from qmworks import (Settings, templates, run, rdkitTools)
from noodles import gather

# User Defined imports
from qmworks.packages.SCM import adf
import plams
# ========== =============

plams.init()

template = "Cc1cc(C)cc(C)c1[PH+](c2c(C)cc(C)cc2C)C[BH3-]"

list_of_modifications = {"Me": "[#0:1]>>[CH3:1]",
                         "Ph": "[#0:1]>>[C:1]c1ccccc1",
                         "tBu": "[#0:1]>>[C:1](C)(C)C",
                         "PhF5": "[#0:1]>>[C:1]c1c(F)c(F)c(F)c(F)c1F",
                         "Mes": "[#0:1]>>[C:1]c1c(C)cc(C)cc1C",
                         "PhCF3_2": "[#0:1]>>[C:1]c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
                         "Cl": "[#0:1]>>[Cl:1]",
                         "CF3": "[#0:1]>>[C:1](F)(F)F"}

HH = plams.Molecule("H_Mes2PCBH2_TS3series1.xyz")
HH.guess_bonds()
newmol = rdkitTools.apply_template(HH, template)
# Change the 2 hydrogens to dummy atoms
temp = rdkitTools.modify_atom(newmol, 47, '*')
temp = rdkitTools.modify_atom(temp, 48, '*')
# Put coordinates of dummy atoms to 0.0, this will force generating new
# coordinates for the new atoms
temp.atoms[47].move_to((0.0, 0.0, 0.0))
temp.atoms[48].move_to((0.0, 0.0, 0.0))

# Replace dummy atoms by new groups, generate coordinates only for new atoms
job_list = []
for mod, smirks in list_of_modifications.items():
    print(mod)
    temp2 = rdkitTools.apply_smirks(temp, smirks)[0]
    temp2 = rdkitTools.apply_smirks(temp2, smirks)[0]
    freeze = rdkitTools.gen_coords(temp2)
    # rdkitTools.write_molblock(temp2)
    # generate job
    s = Settings()
    s.freeze = [a + 1 for a in freeze]
    partial_geometry = adf(templates.geometry.overlay(s), temp2,
                           job_name="partial_opt_" + mod)
    job_list.append(adf(templates.ts, partial_geometry.molecule,
                        job_name="ts_" + mod))

results = run(gather(*job_list), n_processes=1)
