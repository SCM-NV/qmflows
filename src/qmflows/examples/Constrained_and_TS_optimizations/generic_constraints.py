"""An example workflow for generic constrained geometry optimizations."""

__all__ = ['example_generic_constraints']

import scm.plams.interfaces.molecule.rdkit as molkit
from noodles import gather
from qmflows import (dftb, adf, orca, run, Settings, templates, logger)


def example_generic_constraints(*args, **kwargs):
    """Run different job with geometric constrains.

    This examples illustrates that, by using generic keywords, it is possible
    to call different packages interchangeably with the same Settings.
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
            constraint_opt = package(templates.geometry.overlay(
                s), hydrogen_fluoride, job_name)
            jobs.append(constraint_opt)

    # run the jobs
    results = run(gather(*jobs), *args, **kwargs)
    energies = [r.energy for r in results]
    names = [r.job_name for r in results]

    # put resulting energies into a dictionary
    table = {'dftb': {}, 'adf': {}, 'orca': {}}
    for r in results:
        package, distance = r.job_name.split('_')
        table[package][distance] = round(r.energy, 6)

    # print table
    hartree_to_kcalpermol = 627.094
    for package in ['dftb', 'adf', 'orca']:
        row = [package]
        for distance in ['1.0', '1.1', '1.2']:
            val = table[package][distance] - table[package]['1.0']
            row.append(round(val * hartree_to_kcalpermol, 2))
        logger.info('{:10s} {:10.2f} {:10.2f} {:10.2f}'.format(*row))

    return names, energies


if __name__ == "__main__":
    example_generic_constraints()
