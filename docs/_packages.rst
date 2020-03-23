qmflows.packages
================

The :class:`~qmflows.packages.packages.Package` (sub-)classes of QMFlows.

Package-related Functions
-------------------------
.. currentmodule:: qmflows.packages.packages
.. autosummary::
    run

The Package Class
-----------------
.. autosummary::
    Package
    Package.__init__
    Package.__repr__
    Package.__call__
    Package.prerun
    Package.run_job
    Package.postrun
    Package.generic2specific
    Package.get_generic_mapping
    Package.handle_special_keywords

Package Subclasses
------------------
.. currentmodule:: qmflows.packages
.. autosummary::
    ~SCM.ADF
    ~SCM.DFTB
    ~cp2k_package.CP2K
    ~cp2k_mm.CP2KMM
    ~orca.ORCA
    ~package_wrapper.PackageWrapper

Package Instances
------------------
.. autosummary::
    ~SCM.adf
    ~SCM.dftb
    ~cp2k_package.cp2k
    ~cp2k_mm.cp2k_mm
    ~orca.orca

API
---
.. autofunction:: qmflows.packages.packages.run

|

.. autoclass:: qmflows.packages.packages.Package
    :members: generic_mapping, result_type, pkg_name

.. automethod:: qmflows.packages.packages.Package.__init__
.. automethod:: qmflows.packages.packages.Package.__repr__
.. automethod:: qmflows.packages.packages.Package.__call__
.. automethod:: qmflows.packages.packages.Package.prerun
.. automethod:: qmflows.packages.packages.Package.run_job
.. automethod:: qmflows.packages.packages.Package.postrun
.. automethod:: qmflows.packages.packages.Package.generic2specific
.. automethod:: qmflows.packages.packages.Package.get_generic_mapping
.. automethod:: qmflows.packages.packages.Package.handle_special_keywords

|

.. autoclass:: qmflows.packages.SCM.ADF

|

.. autoclass:: qmflows.packages.SCM.DFTB

|

.. autoclass:: qmflows.packages.cp2k_package.CP2K

|

.. autoclass:: qmflows.packages.cp2k_mm.CP2KMM

|

.. autoclass:: qmflows.packages.orca.ORCA

|

.. autoclass:: qmflows.packages.package_wrapper.PackageWrapper
    :noindex:

|

.. autodata:: qmflows.packages.SCM.adf
.. autodata:: qmflows.packages.SCM.dftb
.. autodata:: qmflows.packages.cp2k_package.cp2k
.. autodata:: qmflows.packages.cp2k_mm.cp2k_mm
.. autodata:: qmflows.packages.orca.orca
