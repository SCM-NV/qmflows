"""A module with examples related to frequency calculations."""

__all__ = ['example_freqs']

from noodles import gather
from qmflows import (adf, dftb, orca, run, Settings, templates, find_first_job)
import scm.plams.interfaces.molecule.rdkit as molkit


def is_successful(result):
    """Define the condition for a successful calculation."""
    return result.status not in {"failed", "crashed"}


def example_freqs(*args, n_processes=1, **kwargs):
    """An examples which illustrates the possibility using different packages interchangeably.

    Analytical frequencies are not available for B3LYP in ADF.
    This workflow captures the resulting error and submits the same job to ORCA.
    """
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
    results = run(gather(*jobs), *args, n_processes=n_processes, **kwargs)

    # extrac results
    freqs = [r.frequencies[-3:] for r in results]
    functionals = ['pbe', 'b3lyp', 'blyp']

    # Print the result
    table = ["{:10s}{:10.3f}{:10.3f}{:10.3f}\n".format(fun, *fs)
             for fun, fs in zip(functionals, freqs)]
    print(table)

    return freqs


if __name__ == "__main__":
    example_freqs()
