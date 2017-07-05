
Basic usage tutorial
====================

1. Installation in Unix 
-------------------------

-  | conda installation. Type in your console the following command:

.. code:: bash
	  
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

-  | then add miniconda to your path

.. code:: bash
	  
   miniconda.sh -b -p $HOME/miniconda

-  | create new virtual environment

.. code:: bash
	  
   conda create -q -n qmworks

-  | Install dependecies

.. code:: bash
	  
   conda install --name qmworks -c anaconda hdf5
   conda install --name qmworks -c https://conda.anaconda.org/rdkit rdkit

-  | Start environment

.. code:: bash
	  
   source activate qmworks

-  install **qmworks**

   .. code:: bash

        pip install qmworks --upgrade

   **You are ready to start!**

Starting the environment 
--------------------------

Once *QMWORKS* has been installed the user should run the following
command to initialize the environment:

.. code:: bash

    [user@int1 ~]$ source activate qmworks
    discarding /home/user/anaconda3/bin from PATH
    prepending /home/user/anaconda3/envs/qmworks/bin to PATH
    (qmworks)[user@int1 ~]$ python --version
    Python 3.5.2 :: Anaconda custom (64-bit)

To leave the environment the following command is used

    .. code:: bash

    (qmworks)[user@int1 ~]$ source deactivate
    discarding /home/user/anaconda3/envs/qmworks/bin from PATH

For a more comprenhensive discussion of conda environments see:
`https://conda.io/docs/index.html`

What is QMworks?
=================

QMworks is a python library that enables executing complicated workflows
of interdependent quantum chemical (QM) calculations in python. It aims
at providing a common interface to multiple QM packages, enabling easy
and systematic generation of the calculation inputs, as well as
facilitating automatic analysis of the results. Furthermore it is build
on top of the powerful Noodles framework for executing the calculations
in parallel where possible.

The basics: calling packages
=============================

Currently **QMWORKS** offers an interface with the following simulation
software: \* **SCM (ADF and DTFB)** \* **CP2K** \* **ORCA** \*
**GAMESS-US** \* **DIRAC**

 Please make sure that the packages you want to use in QMworks are
installed and active; in most supercomputer the simulation package are
available using a command like (consult your system administrator):

.. code:: bash

    load module superAwesomeQuantumPackage/3.1421

Also some simulation packages required that you configure a ``scratch``
folder. For instance *Orca* requires a **SCR** folder to be defnied
while *ADF* called it **SCM\_TMPDIR**.

 With ``qmworks`` you can write a python script that simply calls one of
the package objects **adf, dftb, cp2k, orca, gamess** or **dirac**. As
arguments to the call, you need to provide a ``settings`` objects
defining the input of a calculation, a molecular geometry, and,
optionally, a job name that enables you to find back the "raw" data of
the calculation later on.

Let's see how this works:

First we define a molecule, for example by reading one from an xyz file:

.. code:: python

    from plams import Molecule
    acetonitrile = Molecule("files/acetonitrile.xyz")
    print(acetonitrile)


.. parsed-literal::

      Atoms: 
        1         C      2.419290      0.606560      0.000000 
        2         C      1.671470      1.829570      0.000000 
        3         N      1.065290      2.809960      0.000000 
        4         H      2.000000      0.000000      1.000000 
        5         H      2.000000      0.000000     -1.000000 
        6         H      3.600000      0.800000      0.000000 
    


Then we can perform geometry optimization on the molecule by a call to
the dftb package object:

.. code:: python

    from qmworks import dftb, templates, run
    job = dftb(templates.geometry, acetonitrile, job_name="dftb_geometry_optimization")
    print(job)


.. parsed-literal::

    <noodles.interface.decorator.PromisedObject object at 0x7f6c8e5a6d30>


As you can see, "job" is a so-called "promised object". It means it
first needs to be "run" by the Noodles scheduler to return a normal
python object.

.. code:: python

    result = run(job, path="tutorial_results", folder="run_one", cache="tutorial_cache.json")
    print(result)


.. parsed-literal::

    [09:14:04] PLAMS working folder: /home/lars/workspace/qmworks/jupyterNotebooks/tutorial_results/run_one
    ╭─(running jobs)
    │ Running dftb dftb_geometry_optimization...
    ╰[s[1A[50C([38;2;60;180;100m✔[0m)[u─(success)
    <qmworks.packages.SCM.DFTB_Result object at 0x7f6c8e30bcf8>


We can easily retrieve the calculated properties from the DFTB
calculation such as the dipole or the optimized geometry for use in
subsequent calculations.

.. code:: python

    print("Dipole: ", result.dipole)
    print(result.molecule)


.. parsed-literal::

    Dipole:  [1.0864213029, -1.9278296041, -0.0]
      Atoms: 
        1         C      2.366998      0.579794     -0.000000 
        2         C      1.660642      1.834189      0.000000 
        3         N      1.089031      2.847969      0.000000 
        4         H      2.100157      0.010030      0.887206 
        5         H      2.100157      0.010030     -0.887206 
        6         H      3.439065      0.764079     -0.000000 
    


 Settings and templates
=======================

In the above example ``templates.geometry`` was actually a predefined
Settings object. You can define and manipulate Settings in a completely
flexible manner as will be explained in this section. To facilitate
combining different packages in one script, QMworks defines a set of
commonly used generic keywords, which can be combined with package
specific keywords, to provide maximum flexibility.

.. code:: python

    from qmworks import Settings
    s = Settings()
    s.basis = "DZP"
    s.specific.adf.basis.core = "large"
    s.freeze = [1,2,3]
    print(s)


.. parsed-literal::

    basis: 	DZP
    freeze: 	[1, 2, 3]
    specific: 	
             adf: 	
                 basis: 	
                       core: 	large
    


This code snippet illustrates that the ``Settings`` can be specified in
two ways, using generic or specific keywords. Generic keywords represent
input properties that are present in most simulation packages like a
*basis set* while *specific* keywords allow the user to apply specific
keywords for a package that are not in a generic dictionary.

 Expert info: *Settings* are a subclass of python
`dictionaries <https://docs.python.org/3.5/tutorial/datastructures.html#dictionaries>`__
to represent herarchical structures, like

In QMworks/PLAMS multiple settings objects can be combined using the
``overlay`` function.

.. code:: python

    merged_settings = templates.geometry.overlay(s)
    print(merged_settings)


.. parsed-literal::

    basis: 	DZP
    freeze: 	[1, 2, 3]
    specific: 	
             adf: 	
                 basis: 	
                       core: 	large
                       type: 	SZ
                 geometry: 	
                          optim: 	delocal
                 integration: 	
                             accint: 	6.0
                 scf: 	
                     converge: 	1e-06
                     iterations: 	100
                 xc: 	
                    __block_replace: 	True
                    lda: 	
             cp2k: 	
                  force_eval: 	
                             dft: 	
                                 basis_set_file_name: 	
                                 mgrid: 	
                                       cutoff: 	400
                                       ngrids: 	4
                                 potential_file_name: 	
                                 qs: 	
                                    method: 	gpw
                                 scf: 	
                                     OT: 	
                                        N_DIIS: 	7
                                        minimizer: 	DIIS
                                        preconditioner: 	full_single_inverse
                                     eps_scf: 	1e-06
                                     max_scf: 	200
                                     scf_guess: 	atomic
                                 xc: 	
                                    xc_functional: 	pbe
                             subsys: 	
                                    cell: 	
                                         periodic: 	xyz
                  global: 	
                         print_level: 	low
                         project: 	qmworks-cp2k
                         run_type: 	geometry_optimization
                  motion: 	
                         geo_opt: 	
                                 max_iter: 	500
                                 optimizer: 	bfgs
                                 type: 	minimization
             dftb: 	
                  dftb: 	
                       resourcesdir: 	DFTB.org/3ob-3-1
                  task: 	
                       runtype: 	GO
             dirac: 	
             gamess: 	
                    basis: 	
                          gbasis: 	n21
                          ngauss: 	3
                    contrl: 	
                           dfttyp: 	pbe
                           runtyp: 	optimize
                           scftyp: 	rhf
             orca: 	
                  basis: 	
                        basis: 	sto_sz
                  method: 	
                         functional: 	lda
                         method: 	dft
                         runtyp: 	opt
    


The *overlay* method merged the template containing default settings for
geometry optimizations with different packages with the arguments
provided by the user

resulting in:

Note that the generic and specific keywords still exist next to each
other and may not be consistent (e.g. different basis sets are defined
in generic and specific keywords). Upon calling a package with a
Settings object, the generic keywords are first translated into package
specific keywords and combined with the relevant user defined specific
keywords. In this step, the settings defined in generic keywords take
preference. Subsequently, the input file(s) for the given package is/are
generated, based on the keywords after **specific.[package]** based on
the `PLAMS software <https://www.scm.com/doc/plams/index.html>`__.

.. code:: python

    from qmworks import adf
    print(adf.generic2specific(merged_settings))


.. parsed-literal::

    basis: 	DZP
    freeze: 	[1, 2, 3]
    specific: 	
             adf: 	
                 basis: 	
                       core: 	large
                       type: 	DZP
                 constraints: 	
                             atom 2: 	
                             atom 3: 	
                             atom 4: 	
                 geometry: 	
                          optim: 	cartesian
                 integration: 	
                             accint: 	6.0
                 scf: 	
                     converge: 	1e-06
                     iterations: 	100
                 xc: 	
                    __block_replace: 	True
                    lda: 	
             cp2k: 	
                  force_eval: 	
                             dft: 	
                                 basis_set_file_name: 	
                                 mgrid: 	
                                       cutoff: 	400
                                       ngrids: 	4
                                 potential_file_name: 	
                                 qs: 	
                                    method: 	gpw
                                 scf: 	
                                     OT: 	
                                        N_DIIS: 	7
                                        minimizer: 	DIIS
                                        preconditioner: 	full_single_inverse
                                     eps_scf: 	1e-06
                                     max_scf: 	200
                                     scf_guess: 	atomic
                                 xc: 	
                                    xc_functional: 	pbe
                             subsys: 	
                                    cell: 	
                                         periodic: 	xyz
                  global: 	
                         print_level: 	low
                         project: 	qmworks-cp2k
                         run_type: 	geometry_optimization
                  motion: 	
                         geo_opt: 	
                                 max_iter: 	500
                                 optimizer: 	bfgs
                                 type: 	minimization
             dftb: 	
                  dftb: 	
                       resourcesdir: 	DFTB.org/3ob-3-1
                  task: 	
                       runtype: 	GO
             dirac: 	
             gamess: 	
                    basis: 	
                          gbasis: 	n21
                          ngauss: 	3
                    contrl: 	
                           dfttyp: 	pbe
                           runtyp: 	optimize
                           scftyp: 	rhf
             orca: 	
                  basis: 	
                        basis: 	sto_sz
                  method: 	
                         functional: 	lda
                         method: 	dft
                         runtyp: 	opt
    


In the case of adf the above keywords result in the following input file
for ADF package:

.. code:: python

    adf_job = adf(merged_settings, acetonitrile, job_name='adf_acetonitrile')
    result = run(adf_job, path="tutorial_results", 
                 folder="run_two", cache="tutorial_cache.json")
    print(open('tutorial_results/run_two/adf_acetonitrile/adf_acetonitrile.in').read())


.. parsed-literal::

    [09:14:04] PLAMS working folder: /home/lars/workspace/qmworks/jupyterNotebooks/tutorial_results/run_two
    ╭─(running jobs)
    │ Running adf adf_acetonitrile...
    [s[1A[50C([38;2;60;180;100m✔[0m)[u╰─(success)
    atoms
          1         C      2.419290      0.606560      0.000000 
          2         C      1.671470      1.829570      0.000000 
          3         N      1.065290      2.809960      0.000000 
          4         H      2.000000      0.000000      1.000000 
          5         H      2.000000      0.000000     -1.000000 
          6         H      3.600000      0.800000      0.000000 
    end
    
    basis
      core large
      type DZP
    end
    
    constraints
      atom 2
      atom 3
      atom 4
    end
    
    geometry
      optim cartesian
    end
    
    integration
      accint 6.0
    end
    
    scf
      converge 1e-06
      iterations 100
    end
    
    xc
      lda
    end
    
    end input
    


 Combining multiple jobs 
=========================

Multiple jobs can be combined, while calling the run function only once.
The script below combines components outlined above:

.. code:: python

    from plams import Molecule
    from qmworks import dftb, adf, templates, run, Settings
    
    acetonitrile = Molecule("files/acetonitrile.xyz")
    
    dftb_opt = dftb(templates.geometry, acetonitrile, job_name="dftb_opt")
    
    s = Settings()
    s.basis = "DZP"
    s.specific.adf.basis.core = "large"
    adf_single = adf(templates.singlepoint.overlay(s), dftb_opt.molecule, job_name="adf_single")
    
    adf_result = run(adf_single, path="tutorial_results", folder="workflow", cache="tutorial_cache.json")
    print(adf_result.molecule)
    print(adf_result.energy)


.. parsed-literal::

    [09:15:08] PLAMS working folder: /home/lars/workspace/qmworks/jupyterNotebooks/tutorial_results/workflow
    ╭─(running jobs)
    │ Running dftb dftb_opt...
    [s[1A[50C([38;2;60;180;100m✔[0m)[u│ Running adf adf_single...
    [s[1A[50C([38;2;60;180;100m✔[0m)[u╰─(success)
      Atoms: 
        1         C      0.000000      0.000000      0.656511 
        2         C      0.000000      0.000000     -0.783088 
        3         N      0.000000      0.000000     -1.946913 
        4         H     -0.512221     -0.887193      1.022016 
        5         H      1.024442      0.000000      1.022016 
        6         H     -0.512221      0.887193      1.022016 
    
    -1.4094874734528888


In this case the second task adf\_single reads the molecule optimized in
the first job dftb\_opt. Note that dftb\_opt as well as
dftb\_opt.molecule are promised objects. When **run** is applied to the
adf\_single job, noodles builds a graph of dependencies and makes sure
all the calculations required to obtain **adf\_result** are performed.

All data related to the calculations, i.e. input files generated by
QMworks and the resulting output files generated by the QM packages are
stored in folders named after the job\_names, residing inside a results
folder:

.. code:: python

    ls tutorial_results


.. parsed-literal::

    [0m[01;34mrun_one[0m/  [01;34mrun_two[0m/  [01;34mworkflow[0m/


.. code:: python

    ls tutorial_results/workflow


.. parsed-literal::

    [0m[01;34madf_single[0m/  [01;34mdftb_opt[0m/  workflow.log


.. code:: python

    ls tutorial_results/workflow/adf_single


.. parsed-literal::

    adf_single.dill  adf_single.in   [0m[01;32madf_single.run[0m*  logfile  t21.H
    adf_single.err   adf_single.out  adf_single.t21   t21.C    t21.N


