# Version 0.10.2

## New
  * Allow other cp2k executable: ``cp2k.sopt``, ``cp2k.psmp``, etc.


# Version 0.10.1 (09/06/2020)

## Changed
  * Exposed ``InitRestart`` to the main QMFlows ``__init__.py`` file.
  * Exchanged ``plams.init()`` / ``plams.finish()`` for ``qmflows.InitRestart`` in the ``qmflows.run()`` function.
  * Store the ``cache.db`` file in the PLAMS working directory.

# Version 0.10.0 (XX/03/2020)

## Added
  * Introduced the ``CP2KMM`` class for classical forcefield calculations with CP2K: [qmflows/pull/150](https://github.com/SCM-NV/qmflows/pull/150).
  * Introduced the ``PackageWrapper`` class: [qmflows/pull/149](https://github.com/SCM-NV/qmflows/pull/149).
  * Introduced updates and code-style improvements to the ``Package`` and ``Result`` classes: [qmflows/pull/146](https://github.com/SCM-NV/qmflows/pull/146)
  * Added workflow for [GitHub Actions](https://github.com/SCM-NV/qmflows/actions)

## Removed
  * [Removed references to Dirac](https://github.com/SCM-NV/qmflows/issues/152)
  * [Removed Pymonad](https://github.com/SCM-NV/qmflows/issues/156)
  * Remove support for [FDE](https://github.com/SCM-NV/qmflows/issues/171)

# Changed
  * Used [Path](https://github.com/SCM-NV/qmflows/issues/153) instead of ``str``.


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
