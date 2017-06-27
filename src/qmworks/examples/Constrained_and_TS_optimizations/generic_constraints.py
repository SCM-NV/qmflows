__all__ = ['example_generic_constraints']

from noodles import gather
from qmworks import dftb, adf, orca, run, Settings, templates, molkit
import numpy as np

def example_generic_constraints():
    """
    This examples illustrates that, by using generic keywords, it is possible
    to call different packages interchangeably with the same Settings
    """

    # build hydrogen fluoride molecule
    hydrogen_fluoride = molkit.from_smiles('F[H]')

    # loop over distances
    jobs = []
    for distance in [1.0, 1.1, 1.2]:
        s = Settings()
        s.constraint['dist 1 2'] = distance
        # loop over packages
        for package in [dftb, adf, orca]:
            job_name = package.pkg_name + '_' + str(distance)
            constraint_opt = package(templates.geometry.overlay(s), hydrogen_fluoride, job_name)
            jobs.append(constraint_opt)

    # run the jobs
    results = run(gather(*jobs))
    energies = np.array([r.energy for r in results])
    names = [r.job_name for r in results]
    
    return names, energies
