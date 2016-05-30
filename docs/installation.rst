
Installation
============

**QMWorks** depends on a set of libraries that are not pat of the python ecosystem.
First you need to install these dependecies_ then you can follow the installation
of the library using a virtual-environment_ or you can install the python dependecies
from source_.


.. _dependecies:
Dependencies installation
-------------------------

- Install the HDF5 library (note that this library is installed
  by default in most clusters),
 
    - Ubuntu: ``sudo apt-get install libhdf5-dev``

    - Mac (Brew): ``brew install hdf5``

    - Mac (Ports): ``port install hdf5``

    - If you are using linux you should install several libraries before installing python,

      Ubuntu: ``sudo apt-get install build-essential libssl-dev openssl libblas-dev liblapack-dev``
    

.. _virtual-environment:
Installation using a virtual environment (recommended)
------------------------------------------------------

- Download miniconda for python 3.5: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_. 

- Reopen terminal (or type ``source ~/.bashrc``).

- Create a new virtual environment using the following commands:

  - ``conda create -n qmworks python`` 

- Active the new virtual environment
  
  - ``source activate qmworks``

- install rdkit using the following command:

  - ``conda install -y -q -c https://conda.anaconda.org/rdkit rdkit``
    
- Create a new directory *escience* in your home folder.

- Move to the *escience* folder.
  
- Clone the packages using the following commands:
  
   - ``git clone -b escience git@gitlab.pyadf.org:e-science/plams.git``
   - ``git clone git@gitlab.pyadf.org:e-science/workflow-engine.git``    
   - ``git clone git@gitlab.pyadf.org:e-science/qmworks.git``

- Change to directory *$HOME/escience/plams* and type:
  
  ``pip install .``

- Add the following line to the ``$HOME/miniconda3/envs/qmworks/bin/activate`` file,
  
  ``export PLAMSDEFAULTS=$HOME/escience/plams/utils/plams_defaults.py``

- Type in your terminal,

  ``source activate qmworks``
  
- Change to directory *$HOME/escience/workflow-engine* and type:

  ``pip install .``  

- Change to directory *$HOME/escience/qmworks* and type:
  
  ``pip install .``

Now you are ready to use *qmworks*. To exit the virtual environment type  ``source deactivate``.
 

**Notes:**

- Once the libraries and the virtual environment are installed, you only need to type
  ``source active qmworks`` each time that you want to use the software.
 
.. _optional:
Optional packages
-----------------
If you want to use the pyXenon_ interface available with noodles you need to,

- Change to directory *$HOME/escience* and clone the *pyxenon* repository (you must have a github account):

  -``git clone git@github.com:NLeSC/pyxenon.git``

- Change to directory *$HOME/escience/pyxenon* and type:

  -``make install``

- Notice that in order to use xenon you need either a Software Development Kit (**SDK**) from Sun/Oracle or
  Java Development Kit (JDK_). Usually the JDK is install by default in most of the supercomputers. If you
  are using *Xenon* to communicate with remote server from your local machine, see JDK_.
  
.. _source:  
Installation from source
------------------------


- Download and install *python3.5* from: download_.


- Download RDKit from sourceforge_, extract it and follow the instructions inside the RDKit
  folder at *Docs/Book/Install.rst*
    

- Clone the packages using the following commands:
  
   - ``git clone -b escience git@gitlab.pyadf.org:e-science/plams.git``
   - ``git clone git@gitlab.pyadf.org:e-science/workflow-engine.git``    
   - ``git clone git@gitlab.pyadf.org:e-science/qmworks.git``

- Change to directory *$HOME/escience/plams* and type:
  
  ``pip install . --user``

- Add the following line to the ``$HOME/.bash_profile`` (or ``$HOME/.bashrc``) file,
  
  ``export PLAMSDEFAULTS=$HOME/escience/plams/utils/plams_defaults.py``

- Type in your terminal,

  ``source ~/.bash_profile``
  
- Change to directory *$HOME/escience/workflow-engine* and type:

  ``pip install . --user``  

- Change to directory *$HOME/escience/qmworks* and type:
  
  ``pip install . --user``

  
  

.. _miniconda: http://conda.pydata.org/miniconda.html

.. _anaconda: https://www.continuum.io/downloads

.. _installConda: http://conda.pydata.org/docs/install/quick.html

.. _download: https://www.python.org/downloads/

.. _sourceforge: https://sourceforge.net/projects/rdkit/files/rdkit/

.. _JDK: http://www.oracle.com/technetwork/java/javase/downloads/index.html

.. _pyXenon: https://github.com/NLeSC/pyxenon
