__all__ = ['example_freqs']
"""
This examples illustrates the possibility to use different packages interchangeably.
Analytical frequencies are not available for B3LYP in ADF
This workflow captures the resulting error and submits the same job to ORCA.
"""
from noodles import gather
from qmworks import dftb, adf, orca, run, Settings, templates, molkit, find_first_job
import numpy as np

def is_successful(result):
    """
    Define the condition for a successful calculation
    """
    return result.status not in ["failed", "crashed"]


def example_freqs():
    # Generate water molecule
    water = molkit.from_smiles('[OH2]', forcefield='mmff')

    # Pre-optimize the water molecule
    opt_water = dftb(templates.geometry, water, job_name="dftb_geometry")

    jobs = []

    # Generate freq jobs for 3 functionals
    for functional in ['pbe', 'b3lyp', 'blyp']:
        s = Settings()
        s.basis = 'DZ'
        s.functional = functional
        # Try to perform the jobs with adf or orca, take result from  first successful calculation
        freqjob = find_first_job(is_successful, [adf, orca], templates.freq.overlay(s),
                                 opt_water.molecule, job_name=functional)
        jobs.append(freqjob)

    # Run workflow
    results = run(gather(*jobs), n_processes=1)

    freqs = np.array([r.frequencies[-3:] for r in results])

    return freqs

if __name__ == "__main__":
    example_freqs()
