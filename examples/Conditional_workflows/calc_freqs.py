from noodles import gather, schedule, quote, unquote
from qmworks import dftb, adf, orca, run, Settings, templates, molkit
# from qmworks import orca as adf

# This examples illustrates that the possibility to use different packages
# interchangeably. Analytical frequencies are available not for some functionals in ADF
# This workflow captures the resulting error and submits the same job to ORCA

water = molkit.from_smiles('[OH2]')

hartree_to_kcalpermol = 627.094

def is_successful(result):
    return result.status not in ["failed", "crashed"]

def find_first(pred, lst):
    """Receives a predicate (non-scheduled) and a list of promised objects,
    Promises are executed until we find one for which the predicate returns
    true."""
    if lst:
        return s_find_first(pred, lst[0], [quote(l) for l in lst[1:]])
    else:
        return None

@schedule
def s_find_first(pred, first, lst):
    if pred(first):
        return first
    elif lst:
        return s_find_first(pred, unquote(lst[0]), lst[1:])
    else:
        return None



opt_water = dftb(templates.geometry, water, job_name="dftb_geometry")

jobs = []

# loop over functionals
for functional in ['pbe', 'b3lyp', 'blyp']:
    s=Settings()
    s.functional = functional
    adfjob = adf(templates.freq.overlay(s), opt_water.molecule, job_name=functional)
    orcajob = orca(templates.freq.overlay(s), opt_water.molecule, job_name=functional)
    freqjob = find_first(is_successful, [adfjob, orcajob])
    jobs.append(freqjob)

results = run(gather(*jobs), n_processes=1)

for i in range(3):
    print("{:10s}{:10.3f}{:10.3f}{:10.3f}".format(['pbe', 'b3lyp', 'blyp'][i], *results[i].frequencies[-3:]))