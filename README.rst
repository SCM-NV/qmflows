
.. image:: https://github.com/SCM-NV/qmflows/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/SCM-NV/qmflows/actions
.. image:: https://codecov.io/gh/SCM-NV/qmflows/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/SCM-NV/qmflows
.. image:: https://readthedocs.org/projects/qmflows/badge/?version=latest
   :target: https://qmflows.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3274284.svg
   :target: https://doi.org/10.5281/zenodo.3274284
.. image:: https://badge.fury.io/py/qmflows.svg
   :target: https://badge.fury.io/py/qmflows
.. image:: qmflows.png

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
 4. Jobs failure detection and recovery.
 5. Numerical data storage (h5py_).

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


- install rdkit_ and h5py_ using conda:

  - ``conda install -y -q -c conda-forge rdkit h5py``

  - Note that ``rdkit`` is optional for Python 3.7 and later.

.. _installation:

Package installation
--------------------
Finally install the package:

- Install **QMFlows** using pip:
  - ``pip install qmflows``

Now you are ready to use *qmflows*.


  **Notes:**

  - Once the libraries and the virtual environment are installed, you only need to type
    ``conda activate qmflows`` each time that you want to use the software.


.. _documentation: https://qmflows.readthedocs.io/en/latest/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _anaconda: https://www.anaconda.com/distribution/#download-section
.. _installConda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _Noodles: http://nlesc.github.io/noodles/
.. _h5py: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _jupyter-notebook: http://jupyter.org/
.. _tutorial-qmflows: https://github.com/SCM-NV/qmflows/tree/master/jupyterNotebooks
.. _examples: https://github.com/SCM-NV/qmflows/tree/master/src/qmflows/examples
.. _PLAMS: https://github.com/SCM-NV/PLAMS
