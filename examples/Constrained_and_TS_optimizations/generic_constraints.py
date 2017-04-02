from noodles import gather
from qmworks import dftb, adf, orca, run, Settings, templates, molkit

# This examples illustrates that, by using generic keywords, it is possible
# to call different packages interchangeably with the same Settings

# build hydrogen fluoride molecule
hydrogen_fluoride = molkit.from_smiles('F[H]')

hartree_to_kcalpermol = 627.094

jobs = []

# loop over distances
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

# put resulting energies into a dictionary
table = {'dftb':{}, 'adf':{}, 'orca':{}}
for r in results:
    package, distance = r.job_name.split('_')
    table[package][distance] = round(r.energy, 6)

# print table
for package in ['dftb', 'adf', 'orca']:
    row = [package]
    for distance in ['1.0', '1.1', '1.2']:
        row.append(round((table[package][distance] - table[package]['1.0']) * hartree_to_kcalpermol, 2))
    print('{:10s} {:10.2f} {:10.2f} {:10.2f}'.format(*row))


