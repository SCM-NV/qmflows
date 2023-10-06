# Version 1.0.0 (6/10/2023)

## New
 * Add new aliases for CP2K potentials

## Changed
 * Dropped support for Python 3.6 and 3.7
 * Finished the QMFlows 0.12.0 deprecations
 * Switched from setup.cfg to pyproject.toml
 * Misc annotation improvements
 * Switch from setup.py to pyproject.toml
 * Add support for Python 3.12


# Version 0.13.0 (19/04/2023)

## New
 * Revert the `qmflows.packages.registry` deprecation
 * Add formal python 3.11 support
 * General CI clean up

## Changed
 * Changed the default CP2K executable extension from .popt to .ssmp


# Version 0.12.1 (18/05/2022)

## New
 * Add the MO index and occupation numbers to the CP2K orbital output.

## Changed
 * Explicitly raise when the line with the number of orbitals doesn't have any actual orbitals.

## Fixed
 * Fixed various bugs related to the parsing of unrestricted orbitals.


# Version 0.12.0 (13/04/2022)

## New
 * Add support for parsing CP2K basis sets consisting of multiple exponent sets.

## Changed
 * Removed the upper version bounds of ``noodles`` and ``pyparsing``.
 * Remove MO padding when the requested MO range in CP2K is larger than the number of available MOs.
 * Restructured the layout of the QMFlows' submodules; the old layout has been deprecated.
  * The nested sub-modules ``qmflows.packages``, ``.parsers`` and ``.templates``
    have been flatened somewhat.
  * Removal of the ``qmflows.settings`` module; the `Settings` class should
    be imported from the main ``qmflows`` namespace.
  * Renaming the ``qmflows.parsers.cp2k.cp2KBasisParser`` function to ``cp2k_basis_parser``.
  * Removal of the ``registry``, ``load_properties``, ``SerMolecule`` and ``SerSettings``
    functions from the `qmflows.packages` namespace.
  * Moving the ``example_freqs``, ``example_H2O2_TS``, ``example_generic_constraints`` and
    ``example_partial_geometry_opt`` functions from ``qmflows`` to the ``qmflows.examples`` namespace.

## Fix
 * Do not extract the CP2K version via ``cp2k.popt --version``, read it directly from the .out file.


# Version 0.11.2 (07/01/2022)

## Changed
 * Reworked the CP2K basis set parser; allow it to return aliases of basis set names.


# Version 0.11.1 (21/01/2022)

## New
 * Introduced a new template for frequency analyses with CP2K (#278).
 * Allow ``dir()`` to work on result-based generic properties.
 * Added the ``basis`` and ``potential`` generic keywords to CP2K.


# Version 0.11.0 (17/11/2021)

## New
 * Add support for reading CP2K MOs from unrestricted calculations.
 * Add support for reading CP2K >=8.2 MOs.
 * Add a template for (CP2K) cell optimizations: ``qmflows.cell_opt``.
 * Add a generic keyword for the CP2K GAL19 non-bonded forcefield.
 * Add 6 new generic properties to ``qmflows.cp2k`` and ``qmflows.cp2k_mm`` outputs:
   * ``volume``
   * ``forces``
   * ``coordinates``
   * ``temperature``
   * ``lattice``
   * ``pressure``

## Changed
 * Make ``qmflows.Package`` instance more compatible with builtin functions.
 * Remove the unused ``__block_replace`` functionality.
 * Remove the cell parameters from the ``qmflows.cp2k_mm`` templates.
 * Remove the 2-digit restriction from CP2K cell parameters.
 * Check for duplicate keys when parsing .yaml inputs.
 * QMFlows templates are now always copied when getting them (requires Python >= 3.7).
 * Make RDKit an optional dependency (requires Python >= 3.7).

## Fix
 * Fix the ``ResultWrapper`` parameters being ordered incorrectly.
 * Fix ``qmflows.cp2m_mm`` ignoring the ``executable`` key.
 * Fix ``qmflows.InitRestart`` failing on consecutive calls.
 * Fix ``qmflows.CP2KMM_Result`` not inheriting from ``qmflows.CP2K_Result``.
 * Remove usage of the CP2K ``USE_ELEMENT_AS_KIND`` keyword.


# Version 0.10.4 (07/09/2020)
## New
 * Introduced a flag for keeping the Log files

## Fix
 * Improve CP2K error reporting (#209)


# Version 0.10.3 (12/06/2020)

## New
  * Added tests for generating the Sphinx documentation.

## Changed
  * Replaced ``requirements.txt`` with ``.readthedocs.yml``.
  * Fixed the jupyter notebook in the documentation.


# Version 0.10.2 (12/06/2020)

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

 * Used [pyyaml](https://pyyaml.org/wiki/PyYAMLDocumentation) for the [templates](https://github.com/SCM-NV/qmflows/blob/master/src/qmflows/templates/_templates.py) instead of *JSON*
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
