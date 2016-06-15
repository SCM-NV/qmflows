.. packages_:
   
Packages
========

.. currentmodule:: qmworks.packages

The base class |Package| is the library core, it provides the general scaffold to call a quantum code. On top of this infrastructure it has been created several subclasses that contain the specific details for each quantum code. The available interfaces to quantum codes are:

* :class:`~qmworks.packages.SCM.ADF`
* :class:`~qmworks.packages.SCM.DFTB`
* :class:`~qmworks.packages.dirac.DIRAC`
* :class:`~qmworks.packages.cp2k.CP2K`
* :class:`~qmworks.packages.orca.ORCA`



**This class must not be call directly**, instead the correspoding class for the quantum package should be called or in case that there is not an interface to your quantum code, you can make a new subclass that implement the following method:

* ``run_job`` -- This methods takes a |Settings| object a molecule and call a function to create the input automatically and takes cares of the bookkeeping associated with creating new folders, calling the package and retrieving and Object-result.


Instead of implementing all the runners and the bookkeeping functions ourselves, we use the plams_ library. 

For example to carry out a simulation using ADF we call plams as follows ::

  result = plams.ADFJob(molecule=mol, settings=adf_settings).run()
  
See :ref: running for more details about job execution.

 

API
~~~

.. autoclass:: Package
   :members:
   :special-members:
   :exclude-members: __weakref__

SCM
~~~

.. autoclass:: qmworks.packages.SCM.ADF
.. automethod::qmworks.packages.SCM.ADF.molecule
   :members:	      


.. autoclass:: qmworks.packages.SCM.ADF_Result
   :members:	      

.. autoclass:: qmworks.packages.SCM.DFTB
   :members:	      

.. autoclass:: qmworks.packages.SCM.DFTB_Result
   :members:	      

DIRAC
~~~~~
.. autoclass:: qmworks.packages.dirac.DIRAC
   :members:	      

.. autoclass:: qmworks.packages.dirac.DIRAC_Result
   :members:	      
     
CP2k
~~~~

.. autoclass::  qmworks.cp2k_package.CP2K
   :members:

.. autoclass:: qmworks.cp2k_package.CP2K_Result
   :members:

ORCA
~~~~

.. autoclass:: qmworks.packages.orca.ORCA
   :members:





.. _plams: 
