.. packages_:
   
Packages
========

.. currentmodule:: qmflows.packages

The base class |Package| is the library core, it provides the general scaffold to call a quantum code. On top of this infrastructure it has been created several subclasses that contain the specific details for each quantum code. The available interfaces to quantum codes are:

* :class:`~qmflows.packages.SCM.ADF`
* :class:`~qmflows.packages.SCM.DFTB`
* :class:`~qmflows.packages.cp2k_package.CP2K`
* :class:`~qmflows.packages.dirac.DIRAC`
* :class:`~qmflows.packages.gamess.GAMESS`
* :class:`~qmflows.packages.orca.ORCA`



**This class must not be call directly**, instead the correspoding class for the quantum package should be called or in case that there is not an interface to your quantum code, you can make a new subclass that implement the following method:

* ``run_job`` -- This methods takes a |Settings| object a molecule and call a function to create the input automatically and takes cares of the bookkeeping associated with creating new folders, calling the package and retrieving and Object-result.


Instead of implementing all the runners and the bookkeeping functions ourselves, we use the plams_ library. 

For example to carry out a simulation using ADF we call plams as follows ::

  result = plams.ADFJob(molecule=mol, settings=adf_settings).run()
