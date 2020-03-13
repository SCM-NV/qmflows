
# Version 0.10.0 (XX/03/2020)

## Added
  * Introduced the ``CP2KMM`` class for classical forcefield calculations with CP2K: [qmflows/pull/150](https://github.com/SCM-NV/qmflows/pull/150).
  * Introduced the ``PackageWrapper`` class: [qmflows/pull/149](https://github.com/SCM-NV/qmflows/pull/149).
  * Introduced updates and code-style improvements to the ``Package`` and ``Result`` classes: [qmflows/pull/146](https://github.com/SCM-NV/qmflows/pull/146).


# Version 0.9.0 (27/11/2019)

## Changed
  * Use autopep to format the code

## Removed
  * Interface to HDF5
  * Turbomol Parser
  * graphviz dependency


# Version 0.8.0 (17/06/2019)

## Changed

 * Used [pyyaml](https://pyyaml.org/wiki/PyYAMLDocumentation) for the [templates](https://github.com/SCM-NV/qmflows/blob/master/src/qmflows/templates/templates.py) instead of *JSON*
 * Updated documentation
 * Test wiht python 3.7


# Version 0.4.0 (25/02/2019)

## Changed

  * Moved all the functionality to build and analysis quantum dot structures to [their own repo](https://github.com/BvB93/CAT)
  * Moved `molkit` functionality to [PLAMS](https://github.com/SCM-NV/PLAMS)

# 08/01/2019

## Removed
*  Quantum Dots builder functionality moved to [CAT](https://github.com/BvB93/CAT)



# 19/10/2018

## Added
 * Quantum Dots builder functionality

## Changed

 * Use [noodles==0.3.0](https://github.com/NLeSC/noodles/releases)
 * Replace [nose](https://nose.readthedocs.io/en/latest/) with [pytest](https://docs.pytest.org/en/latest/)
 * Imported only the core functionality of the library
 * Used [intra-package-references](https://docs.python.org/3/tutorial/modules.html#intra-package-references)
 * Used `__all__` to limit exposed fuctionality

## Removed

 * Dead code related to all noodles API
 * All the `import *`
 * Dead code from components including PES

## Fixed

 * Job manager issue when removing a SCM job



# 02/12/2018

## Added
 * Ligand MOPAC+COSMO-RS property calculation
 * Inter-ligand activation strain analysis have been (UFF)

## Changed

 * Ligand optimization has been overhauled
