from noodles import gather
from qmworks import dftb, adf, orca, run, Settings, templates, molkit, find_first_job

# This examples illustrates that the possibility to use different packages
# interchangeably. Analytical frequencies are available not for some functionals in ADF
# This workflow captures the resulting error and submits the same job to ORCA

water = molkit.from_smiles('[OH2]')

hartree_to_kcalpermol = 627.094

def is_successful(result):
    return result.status not in ["failed", "crashed"]

opt_water = dftb(templates.geometry, water, job_name="dftb_geometry")

jobs = []

# loop over functionals
for functional in ['pbe', 'b3lyp', 'blyp']:
    s=Settings()
    s.basis = 'DZ'
    s.functional = functional
    freqjob = find_first_job(is_successful, [adf, orca], templates.freq.overlay(s),
                             opt_water.molecule, job_name=functional)
    jobs.append(freqjob)

results = run(gather(*jobs), n_processes=1)

table = ""
for i in range(3):
    table += "{:10s}{:10.3f}{:10.3f}{:10.3f}\n".\
        format(['pbe', 'b3lyp', 'blyp'][i], *results[i].frequencies[-3:])
print(table)