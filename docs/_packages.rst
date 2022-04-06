qmflows.packages
================

The :class:`~qmflows.packages.Package` (sub-)classes of QMFlows.

Package-related Functions
-------------------------
.. currentmodule:: qmflows.packages
.. autosummary::
    run

The Package Class
-----------------
.. autosummary::
    Package
    Package.__init__
    Package.__call__
    Package.prerun
    Package.run_job
    Package.postrun
    Package.generic2specific
    Package.handle_special_keywords

.. autosummary::
    ADF
    DFTB
    CP2K
    CP2KMM
    ORCA
    PackageWrapper

Package Instances
------------------
.. autosummary::
    adf
    dftb
    cp2k
    cp2k_mm
    orca

The Result Class
-----------------
.. autosummary::
    Result
    Result.__init__
    Result.get_property
    Result.results

.. autosummary::
    ADF_Result
    DFTB_Result
    CP2K_Result
    CP2KMM_Result
    ORCA_Result
    ResultWrapper

API
---
.. autofunction:: qmflows.packages.run

|

.. autoclass:: qmflows.packages.Package
    :members: generic_mapping, result_type, pkg_name

.. automethod:: qmflows.packages.Package.__init__
.. automethod:: qmflows.packages.Package.__call__
.. automethod:: qmflows.packages.Package.prerun
.. automethod:: qmflows.packages.Package.run_job
.. automethod:: qmflows.packages.Package.postrun
.. automethod:: qmflows.packages.Package.generic2specific
.. automethod:: qmflows.packages.Package.handle_special_keywords

|

.. autoclass:: qmflows.packages.ADF

|

.. autoclass:: qmflows.packages.DFTB

|

.. autoclass:: qmflows.packages.CP2K

|

.. autoclass:: qmflows.packages.CP2KMM

|

.. autoclass:: qmflows.packages.ORCA

|

.. autoclass:: qmflows.packages.PackageWrapper
    :noindex:

|

.. autoclass:: qmflows.packages.Result
    :members: prop_mapping

.. automethod:: qmflows.packages.Result.__init__
.. automethod:: qmflows.packages.Result.get_property
.. autoproperty:: qmflows.packages.Result.results

|

.. autoclass:: qmflows.packages.ADF_Result
.. autoclass:: qmflows.packages.DFTB_Result
.. autoclass:: qmflows.packages.CP2K_Result
.. autoclass:: qmflows.packages.CP2KMM_Result
.. autoclass:: qmflows.packages.ORCA_Result
.. autoclass:: qmflows.packages.ResultWrapper
    :noindex:

|

.. autofunction:: qmflows.packages.adf
.. autofunction:: qmflows.packages.dftb
.. autofunction:: qmflows.packages.cp2k
.. autofunction:: qmflows.packages.cp2k_mm
.. autofunction:: qmflows.packages.orca
