
.. image:: https://travis-ci.org/SCM-NV/qmflows.svg?branch=master
   :target: https://travis-ci.org/SCM-NV/qmflows 
.. image:: https://img.shields.io/badge/python-3.6-blue.svg

QMFlows
#######
See documentation_ for tutorials and documentation.

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
This library consists of a set of modules written in Python3 to
automate the following tasks:

 1. Input generation.
 2. Handle tasks dependencies (Noodles_).
 3. Advanced molecular manipulation capabilities with (rdkit_).
 4. Numerical data storage and manipulation (HDF5_).
 5. Jobs failure detection and recovery.
 6. Distribution in heterogeneous hardware platforms.    

Tutorial and Examples
---------------------
A tutorial written as a jupyter-notebook_ is available from: tutorial-qmflows_. You can
also access direclty more advanced examples_.
    
 
Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Create a new virtual environment using the following commands:

  - ``conda create -n qmflows`` 

- Activate the new virtual environment
  
  - ``source activate qmflows``

To exit the virtual environment type  ``source deactivate``.
    
    
.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate qmflows``  

Using the conda environment the following packages should be installed:    


- install rdkit_ using conda:

  - ``conda install -y -q --name qmflows -c rdkit rdkit``

- install HDF5_ using conda:

  - ``conda install -y -q --name qmflows -c anaconda h5py``
    
    
.. _installation:

Package installation
--------------------
Finally install the package:

- Update PLAMS_ using pip:
  - ``pip install git+https://github.com/SCM-NV/PLAMS@master#egg=plams-1.2``
    
- Install **QMFlows** using pip:
  - ``pip install git+https://github.com/SCM-NV/qmflows@master#egg=qmflows-0.3.0``

Now you are ready to use *qmflows*.  


  **Notes:**

  - Once the libraries and the virtual environment are installed, you only need to type
    ``conda activate qmflows`` each time that you want to use the software.


.. _Quantum Dots builder:

Quantum Dots builder
--------------------
An example input file (including documentation) is located in QD_input_examples_.
Once the path, input ligands and input cores have been specified the job can be run with ``python qd_input.py``
    
.. _documentation: https://qmflows.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: http://conda.pydata.org/docs/install/quick.html
.. _Noodles: http://nlesc.github.io/noodles/
.. _HDF5: http://www.h5py.org/ 
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _jupyter-notebook: http://jupyter.org/
.. _tutorial-qmflows: https://github.com/SCM-NV/qmflows/tree/master/jupyterNotebooks
.. _examples: https://github.com/SCM-NV/qmflows/tree/master/src/qmflows/examples
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QD_input_examples: https://github.com/SCM-NV/qmflows/blob/master/test/QD_input_examples
