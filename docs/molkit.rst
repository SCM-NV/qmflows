Molkit
======

Molkit is module containing a set of functions to manipulate molecules.

It allows easily building molecules from smiles, for which 3D coordinates are generated automatically.
When specifying a forcefield for geometry refinement, please note that `Ebejer et al. 2012`__
suggested that the MMFF forcefield is better for small molecules with no or a single rotable bond.
For larger molecules UFF appears better.

.. __: http://dx.doi.org/10.1021/ci2004658

Molkit enables modifying existing molecules using reaction smarts.

There is also a function to partition a protein into capped fragments for MFCC_ calculations.

.. _MFCC: http://dx.doi.org/10.1063/1.2906128

API
---

.. automodule:: qmflows.molkit
    :members:
