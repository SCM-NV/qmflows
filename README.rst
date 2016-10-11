.. image:: https://img.shields.io/github/license/SCM-NV/qmworks.svg?maxAge=2592000
   :target: https://github.com/SCM-NV/qmworks/blob/master/LICENSE.md
.. image:: https://travis-ci.org/SCM-NV/qmworks.svg?branch=master
   :target: https://travis-ci.org/SCM-NV/qmworks 
.. image:: https://img.shields.io/badge/python-3.5-blue.svg
.. image:: https://img.shields.io/codacy/grade/e27821fb6289410b8f58338c7e0bc686.svg?maxAge=2592000
   :target: https://www.codacy.com/app/tifonzafel/qmworks/dashboard	
	   
================
QMWorks
================


Motivation
==========
Research on modern computational quantum chemistry relies on a set of computational
tools to carry out calculations. The complexity of the calculations usually requires 
intercommunication between the aforementioned tools, such communication is usually done 
through shell scripts that try to automate input/output actions like: launching 
the computations in a cluster, reading the resulting output and feeding the relevant
numerical result to another program. Such scripts are difficult to maintain and extend,
requiring a significant programming expertise to work with them. Being then desirable a
set of automatic and extensible tools that allows to perform complex simulations in
heterogeneous hardware platforms.

This library tackles the construction and efficient execution of computational chemistry workflows.
This allows computational chemists to use the emerging massively parallel compute environments in
an easy manner and focus on interpretation of scientific data rather than on tedious job submission
procedures and manual data processing. 

Description
===========
This library consists of a set of modules written in Python 3.5 to
automate the following tasks:
 1. Input generation.
 2. Handle tasks dependencies (Noodles_).
 3. Distribution in heterogeneous hardware platforms.
 4. Advanced molecular manipulation capabilities with (rdkit_).
 5. Numerical data storage and manipulation (HDF5_).
 6. Jobs failure detection and recovery.

 
Installation
============
First check that you have available in your system **Python 3.5**. Otherwise, you can download it from here_.
Because **QMWorks** depends on a set of libraries that are not pat of the python ecosystem, you need first
to install these dependecies_  using a virtual-environment_. Finally, you can proceed to the package installation_



.. _virtual-environment:
Installation using a virtual environment (recommended)
======================================================

- Download miniconda for python 3.5: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Reopen terminal (or type ``source ~/.bashrc``).

- Create a new virtual environment using the following commands:

  - ``conda create -n qmworks python=3.5`` 

- Activate the new virtual environment
  
  - ``source activate qmworks``

To exit the virtual environment type  ``source deactivate``.
    
    
.. _dependecies:
Dependencies installation
-------------------------

Using the conda environment the following packages should be installed:    


- install rdkit_ using the following command:

  - ``conda install -y -q --name qmworks -c https://conda.anaconda.org/rdkit rdkit``

- install HDF5_ using conda:

  - ``conda install -y -q --name qmworks -c anaconda hdf5``
    

.. _installation:
Package installation
--------------------
    
- Type in your terminal,

  ``source activate qmworks``  

- Then

  ``pip install https://github.com/SCM-NV/qmworks/tarball/master#egg=qmworks  https://github.com/SCM-NV/plams/tarball/master#egg=plams --upgrade``
  
Now you are ready to use *qmworks*. 
 

**Notes:**

- Once the libraries and the virtual environment are installed, you only need to type
  ``source activate qmworks`` each time that you want to use the software.




.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: http://conda.pydata.org/docs/install/quick.html
.. _Noodles: https://gitlab.pyadf.org/e-science/workflow-engine
.. _HDF5: http://www.h5py.org/ 
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _Plams: https://www.scm.com/documentation/Tutorials/Scripting/first_steps_with_plams/
