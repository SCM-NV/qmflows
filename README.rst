
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
This library consists of a set of modules written in Python 3.5,
offering the following features:

 1. Quantum chemistry Workflows execution using Noodles_.

 2. Automatic input generation for several quantum packages.
 
 3. Broadcasting numerical results between quantum chemistry.

 
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

  - ``conda create -n qmworks python`` 

- Activate the new virtual environment
  
  - ``source activate qmworks``

    
.. _dependecies:
Dependencies installation
-------------------------

Using the conda environment the following packages should be installed:    


- install rdkit using the following command:

  - ``conda install -y -q -c https://conda.anaconda.org/rdkit rdkit``

- install HDF5_ using conda:

  -``conda install -c anaconda hdf5``
    

.. _installation:
Package installation
--------------------

    
- Create a new directory *escience* in your home folder.

- Move to the *escience* folder.
  
- Clone the packages using the following commands:
  
   - ``git clone -b escience git@gitlab.pyadf.org:e-science/plams.git``

- Modify the anaconda activation file for Plams adding the following line ``export PLAMSDEFAULTS=$HOME/escience/plams/utils/plams_defaults.py`` to the file 
  
  ``$HOME/miniconda3/envs/qmworks/bin/activate`` or
  ``$HOME/anaconda/envs/qmworks/bin/activate``
  
  You can find the path to your activation file running the command ``conda info --envs``.

- Type in your terminal,

  ``source activate qmworks``

- Change to directory *$HOME/escience/plams* and type:
  
  ``pip install .``

- Change to directory *$HOME/escience* and type:
  
  -``pip install https://github.com/NLeSC/noodles/tarball/devel#egg=Noodles``
  - ``pip install https://github.com/SCM-NV/qmworks/tarball/master#egg=QuantumWorkflows-0.1.2``

Now you are ready to use *qmworks*. To exit the virtual environment type  ``source deactivate``.
 

**Notes:**

- Once the libraries and the virtual environment are installed, you only need to type
  ``source active qmworks`` each time that you want to use the software.




.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: http://conda.pydata.org/docs/install/quick.html
.. _Noodles: https://gitlab.pyadf.org/e-science/workflow-engine
.. _HDF5: http://www.h5py.org/ 
.. _here: https://www.python.org/downloads/

