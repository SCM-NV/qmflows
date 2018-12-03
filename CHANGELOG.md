
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
