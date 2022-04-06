.. packages_:

Packages
========

The base class :class:`~qmflows.packages.Package` is the library core, it provides the general scaffold to call a quantum code.
On top of this infrastructure it has been created several subclasses that contain the specific details for each quantum code.
The available interfaces to quantum codes are:

* :class:`~qmflows.packages.ADF`
* :class:`~qmflows.packages.DFTB`
* :class:`~qmflows.packages.CP2K`
* :class:`~qmflows.packages.CP2KMM`
* :class:`~qmflows.packages.ORCA`



**This class must not be call directly**, instead the correspoding class for the quantum package should be called or in case that there is not an interface to your quantum code,
you can make a new subclass that implement the following method:

* ``run_job`` -- This methods takes a :class:`~qmflows.settings.Settings` object a molecule and call a function to create the input automatically and takes cares of the bookkeeping associated with creating new folders, calling the package and retrieving and Object-result.


Instead of implementing all the runners and the bookkeeping functions ourselves, we use the plams_ library.

For example to carry out a simulation using ADF we call plams as follows:

.. code:: python

    >>> from scm import plams

    >>> mol = plams.Molecule(...)
    >>> adf_settings = plams.Settings(...)

    >>> result = plams.ADFJob(molecule=mol, settings=adf_settings).run()


.. _plams: https://www.scm.com/doc/plams/index.html
