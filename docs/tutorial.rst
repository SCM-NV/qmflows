
 Table of content 
------------------

1. **Installation**
2. **Packages**
3. **Settings**
4. **Templates**
5. **Molecules**
6. \*\* Runinng a quantum mechanics simulation\*\*
7. **How the run function works?**
8. \*\* Advanced examples\*\*
9. **Exception handling**

 1. Installation in Unix 
-------------------------

-  | conda installation. Type in your console the following command:
   | ``bash    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh``

-  | then add miniconda to your path
   | ``bash    bash miniconda.sh -b -p $HOME/miniconda``

-  create new virtual environment ``bash    conda create -q -n qmflows``

-  Install dependecies
   ``bash     conda install --name qmflows -c anaconda hdf5    conda install --name qmflows -c https://conda.anaconda.org/rdkit rdkit``

-  Start environment ``bash    source activate qmflows``

-  install **qmflows**

   .. code:: bash

        pip install qmflows --upgrade

   \*\* You are ready to start! \*\*

 Starting the environment 
--------------------------

Once *QMFLOWS* has been installed the user should run the following
command to initialize the environment:

.. code:: bash

    [user@int1 ~]$ source activate qmflows
    discarding /home/user/anaconda3/bin from PATH
    prepending /home/user/anaconda3/envs/qmflows/bin to PATH
    (qmflows)[user@int1 ~]$ python --version
    Python 3.5.2 :: Anaconda custom (64-bit)

To leave the environment the following command is used

.. code:: bash

    (qmflows)[user@int1 ~]$ source deactivate
    discarding /home/user/anaconda3/envs/qmflows/bin from PATH

 2. What is QMflows? 
---------------------

QMflows is a python library that enables executing complex workflows of
interdependent quantum chemical (QM) calculations in python. It aims at
providing a common interface to multiple QM packages, enabling easy and
systematic generation of the calculation inputs, as well as facilitating
automatic analysis of the results. Furthermore it is build on top of the
powerful `Noodles <http://nlesc.github.io/noodles/>`__ framework for
executing the calculations in parallel where possible. Together with the
`Plams <https://github.com/SCM-NV/PLAMS>`__ library to interface with
the different QM simulation packages.

 Packages 
----------

Currently ``QMflows`` offers an interface with the following simulation
softwares: \* `SCM <https://www.scm.com/>`__ (ADF and DTFB) \*
`CP2K <https://www.cp2k.org/>`__ \*
`ORCA <https://orcaforum.cec.mpg.de/>`__ \*
`GAMESS-US <http://www.msg.ameslab.gov/gamess/>`__ \*
`DIRAC <http://diracprogram.org/doku.php>`__

If you are interested in having support for other packages, request it
using the `github-issues <https://github.com/SCM-NV/qmflows/issues>`__
system (Sorry but Gaussian is out of the menu!).

With ``qmflows`` you can write a python script that simply calls one of
the package objects **adf, dftb, cp2k, orca, gamess** or **dirac**. As
arguments to the call, you need to provide a ``Settings`` objects
defining the input of a calculation, a molecular geometry represented by
an object called ``Molecule``, and, optionally, a job name that enables
you to find back the "raw" data of the calculation later on.

 Technical note: 
~~~~~~~~~~~~~~~~~

``It is the user responsability to install or load the simulation packages that he/she wants to use. In most supercomputers these simulation packages are available using a command like (consult your system administrator):``

.. code:: bash

    load module superAwesomeQuantumPackage/3.141592

``Also some simulation packages required that you configure a ``scratch`` folder. For instance *Orca* requires a ``SCR`` folder to be defined while *ADF*  called it ``SCM_TMPDIR``.``

 3. QMflows Settings 
---------------------

*Settings* are a subclass of python
`dictionaries <https://docs.python.org/3.5/tutorial/datastructures.html#dictionaries>`__
to represent hierarchical structures, like



.. code:: ipython3

    from qmflows import Settings
    
    s = Settings()
    s.b.z
    s.c.f
    s.c.g = 0

These hierachical resemble the input structure used in most quantum
simulation package. For instance the basis set section in ADF is given
by something like:

::

    Basis
      Type DZP
      Core Large
    End

We can resemble this structure using **Settings**,

.. code:: ipython3

    s = Settings()
    s.specific.adf.basis.basis = "DZP"
    s.specific.adf.basis.core = "Large"

We are creating the *adf* hierarchy under a key called *specific*, this
key is used to differentiate keywords that are unique to a certain
quantum package from those that can be used in several packages as we
will see in the next section.

similarly, we can define ``Settings`` for all the sections

::

    Basis
      Type DZP
    End

    Constraints
      Dist 1 2 1.0
    End

    Geometry
      Optim delocal
    End

    Integration
      Accint 6.0
    End

    Scf
      Converge 1e-06
      Iterations 100
    End

    Xc
      Lda
    End

Represented by the following code

.. code:: ipython3

    s = Settings()
    
    # Basis
    s.specific.adf.basis.basis = "DZP"
    s.specific.adf.basis.core = "Large"
    
    # Constrains
    s.specific.adf.constraints.dist  = "1 2 1.0"
    
    #Geometry
    s.specific.adf.geometry.optim = 'delocal'
    
    #Integration
    s.specific.adf.integration.accint = 6.0
    
    # SCF
    s.specific.adf.scf.converge = 1e-6
    s.specific.adf.scf.iterations = 100
    
    # Functional
    s.specific.adf.xc.lda
    
    print(s)


.. parsed-literal::

    specific: 	
             adf: 	
                 basis: 	
                       basis: 	DZP
                       core: 	Large
                 constraints: 	
                             dist: 	1 2 1.0
                 geometry: 	
                          optim: 	delocal
                 integration: 	
                             accint: 	6.0
                 scf: 	
                     converge: 	1e-06
                     iterations: 	100
                 xc: 	
                    lda: 	
    


You don't need to explicitly declare the ``end`` keyword, *qmflows*
knows how to hande them.

 Generic Keywords 
~~~~~~~~~~~~~~~~~~

Quantum chemistry packages use gaussian type orbitals (GTO) or slater
type orbitals (STO) to perform the simulation. The packages use the same
standards for the basis set and it will be really handy if we can
defined a "generic" keyword for basis sets. Fortunately ``qmflows``
already offers such keyword that can be used among the packages that use
the same basis standard,

.. code:: ipython3

    s = Settings()
    s.basis = "DZP"

Internally **QMflows** will create a hierarchical structure representing
basis *DZP* for the packages that can handle that basis set. Other
generic keyowrds like: ``functional``, ``inithess``, etc. have been
implemented.

The Following table describes some of the available generic keywords and
the Packages where the funcionality is implemented

+------------------+------------------+----------+
| Generic Keyword  | Packages         | Descript |
|                  |                  | ion      |
+==================+==================+==========+
| basis            | ADF, CP2K, Orca  | Set the  |
|                  |                  | Basis    |
|                  |                  | set      |
+------------------+------------------+----------+
| cell\_angles     | CP2K             | Specifie |
|                  |                  | d        |
|                  |                  | the      |
|                  |                  | angles   |
|                  |                  | of the   |
|                  |                  | unit     |
|                  |                  | cell     |
+------------------+------------------+----------+
| cell\_parameters | CP2K             | Specifie |
|                  |                  | d        |
|                  |                  | the      |
|                  |                  | vectors  |
|                  |                  | of the   |
|                  |                  | unit     |
|                  |                  | cell     |
+------------------+------------------+----------+
| constraint       | ADF, Orca        | Constrai |
|                  |                  | n        |
|                  |                  | the      |
|                  |                  | distance |
|                  |                  | ,        |
|                  |                  | angle or |
|                  |                  | dihedral |
|                  |                  | angle    |
|                  |                  | for a    |
|                  |                  | set of   |
|                  |                  | molecule |
|                  |                  | s        |
+------------------+------------------+----------+
| freeze           | ADF, Gamess,     | Freeze a |
|                  | Orca             | set of   |
|                  |                  | atoms    |
|                  |                  | indicate |
|                  |                  | d        |
|                  |                  | by their |
|                  |                  | indexes  |
|                  |                  | or       |
|                  |                  | symbols  |
+------------------+------------------+----------+
| functional       | ADF, CP2K        | Set the  |
|                  |                  | DFT      |
|                  |                  | function |
|                  |                  | al       |
+------------------+------------------+----------+
| inithess         | ADF, Orca        | Provide  |
|                  |                  | an       |
|                  |                  | initial  |
|                  |                  | Hessian  |
|                  |                  | matrix   |
+------------------+------------------+----------+
| optimize         | ADF, DFTB, Orca  | Perform  |
|                  |                  | a        |
|                  |                  | molecula |
|                  |                  | r        |
|                  |                  | optimiza |
|                  |                  | tion     |
+------------------+------------------+----------+
| selected\_atoms  | ADF, Gamess,     | Optimize |
|                  | Orca             | the      |
|                  |                  | given    |
|                  |                  | set of   |
|                  |                  | atoms    |
|                  |                  | while    |
|                  |                  | keeping  |
|                  |                  | the rest |
|                  |                  | fixed    |
+------------------+------------------+----------+
| ts               | ADF, DFTB, Orca  | Carry    |
|                  |                  | out a    |
|                  |                  | transiti |
|                  |                  | on       |
|                  |                  | state    |
|                  |                  | optimiza |
|                  |                  | tion     |
+------------------+------------------+----------+

Note: **QMflows** Does not have chemical intuition and if you provide a
meaningless keyword, like a wrong basis set it will not warm you.

 4. Templates 
--------------

As has been shown so far, **Settings** can be specified in two ways:
generic or specific. Generic keywords represent input properties that
are present in most simulation packages like a *basis set* while
*specific* keywords resemble the input structure of a given package.

*Generic* and *Specific* **Settings** can express both simple and
complex simulation inputs, but it would be nice if we can pre-defined a
set of templates for the most common quantum chemistry simulations like:
single point calculations, geometry optimizations, transition state
optimization, frequency calculations, etc. *qmflows* already has a
pre-defined set of templates containing some defaults that the user can
modify for her/his own purpose. ``Templates`` are stored inside the
``qmflows.templates`` module and are load from *JSON* files. A JSON file
is basically a nested dictionary that is translated to a ``Settings``
object by *qmflows*.

Below it is shown the defaults for single point calculation

.. code:: ipython3

    from qmflows import templates
    templates.singlepoint




.. parsed-literal::

    _ipython_canary_method_should_not_exist_: 	
    specific: 	
             adf: 	
                 basis: 	
                       type: 	SZ
                 integration: 	
                             accint: 	4.0
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
                                 print: 	
                                       mo: 	
                                          add_last: 	numeric
                                          each: 	
                                               qs_scf: 	0
                                          eigenvalues: 	
                                          eigenvectors: 	
                                          filename: 	./mo.data
                                          ndigits: 	36
                                          occupation_numbers: 	
                                 qs: 	
                                    method: 	gpw
                                 scf: 	
                                     added_mos: 	
                                     eps_scf: 	1e-06
                                     max_scf: 	200
                                     scf_guess: 	restart
                                 xc: 	
                                    xc_functional: 	pbe
                             subsys: 	
                                    cell: 	
                                         periodic: 	xyz
                  global: 	
                         print_level: 	low
                         project: 	qmflows-cp2k
                         run_type: 	energy_force
             dftb: 	
                  dftb: 	
                       resourcesdir: 	DFTB.org/3ob-3-1
                  task: 	
                       runtype: 	SP
             dirac: 	
                   DIRAC: 	WAVEFUNCTION
                   HAMILTONIAN: 	LEVY-LEBLOND
                   WAVE FUNCTION: 	SCF
             gamess: 	
                    basis: 	
                          gbasis: 	sto
                          ngauss: 	3
                    contrl: 	
                           dfttyp: 	pbe
                           scftyp: 	rhf
             orca: 	
                  basis: 	
                        basis: 	sto_sz
                  method: 	
                         functional: 	lda
                         method: 	dft




The question is then, *how I can modify a template with my own changes?*

Suppose you are perfoming a bunch of constrained *DFT* optimizations
using ``ADF`` . You need first to define a basis set and the constrains.

.. code:: ipython3

    s = Settings()
    # Basis
    s.basis = "DZP"
    s.specific.adf.basis.core = "Large"
    
    # Constrain
    s.freeze = [1, 2, 3]

We use two *generic* keywords: ``freeze`` to indicate a constrain and
``basis`` to provide the basis set. Also, we introduce an specific
``ADF`` keywords ``core = Large``. Now you merge your **Settings** with
the correspoding template to carry out molecular geometry optimizations,
using a method called ``overlay``.

.. code:: ipython3

    from qmflows import templates
    inp = templates.geometry.overlay(s)

The ``overlay`` method takes as input a template containing a default
set for different packages and also takes the arguments provided by the
user, as shown schematically

This ``overlay`` method merged the defaults for a given packages (*ADF*
in this case) with the input supplied by the user, always given
preference to the user input

Below it is shown a combination of templates, generic and specific
keywords to generate the input for a ``CP2K`` job

.. code:: ipython3

    from qmflows import templates
    
    # Template
    s = templates.singlepoint
    
    # Generic keywords
    s.cell_angles = [90.0, 90.0, 90.0]
    s.cell_parameters=  38.0  
    s.basis = 'DZVP-MOLOPT-SR-GTH'
    s.potential ='GTH-PBE'
    
    # Specific Keywords
    s.specific.cp2k.force_eval.dft.scf.max_scf  = 100
    s.specific.cp2k.force_eval.subsys.cell.periodic = 'None'
    
    print(s)


.. parsed-literal::

    _ipython_canary_method_should_not_exist_: 	
    basis: 	DZVP-MOLOPT-SR-GTH
    cell_angles: 	[90.0, 90.0, 90.0]
    cell_parameters: 	38.0
    potential: 	GTH-PBE
    specific: 	
             adf: 	
                 basis: 	
                       type: 	SZ
                 integration: 	
                             accint: 	4.0
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
                                 print: 	
                                       mo: 	
                                          add_last: 	numeric
                                          each: 	
                                               qs_scf: 	0
                                          eigenvalues: 	
                                          eigenvectors: 	
                                          filename: 	./mo.data
                                          ndigits: 	36
                                          occupation_numbers: 	
                                 qs: 	
                                    method: 	gpw
                                 scf: 	
                                     added_mos: 	
                                     eps_scf: 	1e-06
                                     max_scf: 	100
                                     scf_guess: 	restart
                                 xc: 	
                                    xc_functional: 	pbe
                             subsys: 	
                                    cell: 	
                                         periodic: 	None
                  global: 	
                         print_level: 	low
                         project: 	qmflows-cp2k
                         run_type: 	energy_force
             dftb: 	
                  dftb: 	
                       resourcesdir: 	DFTB.org/3ob-3-1
                  task: 	
                       runtype: 	SP
             dirac: 	
                   DIRAC: 	WAVEFUNCTION
                   HAMILTONIAN: 	LEVY-LEBLOND
                   WAVE FUNCTION: 	SCF
             gamess: 	
                    basis: 	
                          gbasis: 	sto
                          ngauss: 	3
                    contrl: 	
                           dfttyp: 	pbe
                           scftyp: 	rhf
             orca: 	
                  basis: 	
                        basis: 	sto_sz
                  method: 	
                         functional: 	lda
                         method: 	dft
    


 5. Molecule 
-------------

The next component to carry out a simulation is a molecular geometry.
*qmflows* offers a convinient way to read Molecular geometries using the
`Plams <https://www.scm.com/doc/plams/molecule.html>`__ library in
several formats like: *xyz* , *pdb*, *mol*, etc.

.. code:: ipython3

    from plams import Molecule
    acetonitrile = Molecule("files/acetonitrile.xyz")
    print(acetonitrile)

You can also create the molecule one atom at a time

.. code:: ipython3

    from plams import (Atom, Molecule)
    m  = Molecule()
    m.add_atom(Atom(symbol='C', coords=(2.41929, 0.60656 , 0.0)))
    m.add_atom(Atom(symbol='C', coords=(1.67147,  1.82957, 0.0)))
    m.add_atom(Atom(symbol='N', coords=(1.06529, 2.80996, 0.0)))
    m.add_atom(Atom(symbol='H',  coords=(2.0, 0.0, 1.0)))
    m.add_atom(Atom(symbol='H',  coords=(2.0, 0.0, -1.0)))
    m.add_atom(Atom(symbol='H',  coords=(3.6, 0.8, 0.0)))
    print(m)

**QMflows** Can also handle smiles as shown below

.. code:: ipython3

    from qmflows.molkit import from_smiles
    
    # String representing the smile
    smile = 'C1CC2(CCCCC2)C=C1'
    #Molecule creation
    mol = from_smiles(smile)
    print(mol)


.. parsed-literal::

      Atoms: 
        1         C      2.798599     -0.150460      0.089927 
        2         C      1.615862     -0.067227     -0.832063 
        3         C      0.376333      0.019368      0.053118 
        4         C     -0.347606      1.253513     -0.326955 
        5         C     -1.822801      1.252517     -0.204840 
        6         C     -2.446980     -0.058315      0.156076 
        7         C     -1.752081     -1.139264     -0.623091 
        8         C     -0.361478     -1.268307     -0.080395 
        9         C      0.939434      0.095296      1.441284 
       10         C      2.254053      0.292268      1.391042 
       11         H      3.270712     -1.144983      0.118276 
       12         H      3.579260      0.571882     -0.225622 
       13         H      1.583048     -0.956205     -1.489932 
       14         H      1.723819      0.828427     -1.481830 
       15         H     -0.086390      1.553756     -1.374861 
       16         H      0.059603      2.087392      0.301854 
       17         H     -2.168065      1.989676      0.564557 
       18         H     -2.314662      1.598505     -1.150678 
       19         H     -3.508879     -0.050718     -0.200863 
       20         H     -2.505664     -0.256055      1.238793 
       21         H     -2.273284     -2.097029     -0.440451 
       22         H     -1.722856     -0.837576     -1.680683 
       23         H      0.171167     -1.989673     -0.753019 
       24         H     -0.417079     -1.810866      0.897936 
       25         H      0.453939     -0.226576      2.367445 
       26         H      2.901995      0.510654      2.244971 
      Bonds: 
       (1)--1.0--(2)
       (2)--1.0--(3)
       (3)--1.0--(4)
       (4)--1.0--(5)
       (5)--1.0--(6)
       (6)--1.0--(7)
       (7)--1.0--(8)
       (3)--1.0--(9)
       (9)--2.0--(10)
       (10)--1.0--(1)
       (8)--1.0--(3)
       (1)--1.0--(11)
       (1)--1.0--(12)
       (2)--1.0--(13)
       (2)--1.0--(14)
       (4)--1.0--(15)
       (4)--1.0--(16)
       (5)--1.0--(17)
       (5)--1.0--(18)
       (6)--1.0--(19)
       (6)--1.0--(20)
       (7)--1.0--(21)
       (7)--1.0--(22)
       (8)--1.0--(23)
       (8)--1.0--(24)
       (9)--1.0--(25)
       (10)--1.0--(26)
    


The Molecule class has an extensive functionally to carry out molecular
manipulations, for a comprenhesive disccusion about it have a look at
the `molecule
documentation <https://www.scm.com/doc/plams/molecule.html>`__. Also the
module ``qmflows.molkit`` contains an extensive functionality to apply
transformation over a molecule using the
`RDKit <http://www.rdkit.org/>`__ library.

 6. Runinng a quantum mechanics simulation 
-------------------------------------------

We now have our components to perform a calculation: **Settings** and
**Molecule**. We can now invoke a quantum chemistry package to perform
the computation,

.. code:: ipython3

    from qmflows import adf
    optmized_mol_adf = adf(inp, acetonitrile, job_name='acetonitrile_opt')

the previous code snippet *does not execute the code immediatly*,
instead the simulation is started when the user invokes the run
function, as shown below

.. code:: python

    from plams import Molecule
    from qmflows import (adf, run, Settings, templates)

    # Settings
    s = templates.geometry
    s.basis = "DZP"
    s.specific.adf.basis.core = "Large"
    s.freeze = [1, 2, 3]

    # molecule 
    from plams import Molecule
    acetonitrile = Molecule("acetonitrile.xyz")

    # Job 
    optimized_mol_adf = adf(s, acetonitrile, job_name='acetonitrile_opt')

    # run the  job
    result = run(optimized_mol_adf.molecule, folder='acetonitrile')

you can run the previous script by saving it in a file called
``acetonitrile_opt.py`` and typing the following command in your
console:

.. code:: bash

    (qmflows)[user@int1 ~]$ python acetonitrile_opt.py

you will then see in your ``current work directory`` something similar
to the following

.. code:: bash

    (qmflows)[user@int1 ~]$ ls 
    acetonitrile      acetonitrile_opt.py   cache.json   acetonitrile.xyz  

 acetonitrile is the folder containing the output from the quantum
package call (``ADF`` in this case). The ``cache.json`` file contains
all the information required to perform a restart, as we will explore
below. Inside the acetonitrile you can find the input/output files
resulting from the simulation

``bash (qmflows)[user@int1 ~]$ ls acetonitrile   acetonitrile.log  acetonitrile_opt``

``bash (qmflows)[user@int1 ~]$ ls acetonitrile/acetonitrile_opt  acetonitrile_opt.dill  acetonitrile_opt.out  logfile  t21.N  acetonitrile_opt.err   acetonitrile_opt.run  t21.C  acetonitrile_opt.in    acetonitrile_opt.t21  t21.H``

 Extracting Properties 
-----------------------

In general, properties are extracted using the standard
``Object.attribute`` notation in python, as shown below.

.. code:: python

    result = optmized_mol_adf.molecule

Some of the available properties are shown in the following table,

+------------+---------------+-----------------------------------------+
| Property   | type          | Description                             |
+============+===============+=========================================+
| dipole     | Double        | Molecular dipole mopment                |
+------------+---------------+-----------------------------------------+
| energy     | Double        | Total energy                            |
+------------+---------------+-----------------------------------------+
| enthalpy   | Double        | Enthalpy                                |
+------------+---------------+-----------------------------------------+
| gradient   | Numpy array   | First derivatives of the energy         |
+------------+---------------+-----------------------------------------+
| hessian    | Numpy array   | Second derivative of the energy         |
+------------+---------------+-----------------------------------------+
| molecule   | Molecule      | Object representing a physical entity   |
+------------+---------------+-----------------------------------------+
| runtime    | Double        | Time spent in the simulation            |
+------------+---------------+-----------------------------------------+

On the background *QMflows* has a mechanism to read the properties from
the output files and make them available inside Python.

 Communicating different packages 
----------------------------------

We can use the previous optimized geometry for further calculations
using for instance another package like *Orca* to run a frequencies
calculation,

.. code:: ipython3

    from qmflows import orca
    s2 = Settings()
    s2.specific.orca.main = "freq"
    s2.specific.orca.basis.basis = 'sto_sz'
    s2.specific.orca.method.functional = 'lda'
    s2.specific.orca.method.method = 'dft'
    
    job_freq = orca(s2, optmized_mol_adf)
    
    frequencies = job_freq.frequencies

The whole script is

.. code:: python

    from qmflows import (adf, orca, run, templates, Settings)
    from plams import Molecule
    import plams

    def main():
        s = templates.geometry
        s.basis = "DZP"
        s.specific.adf.basis.core = "large"

        acetonitrile = Molecule("files/acetonitrile.xyz")
        job = adf(inp, acetonitrile)
        optmized_mol_adf = job.molecule

        s2 = Settings()
        s2.specific.orca.main = "freq"
        s2.specific.orca.basis.basis = 'sto_sz'
        s2.specific.orca.method.functional = 'lda'
        s2.specific.orca.method.method = 'dft'

        job_freq = orca(s2, optmized_mol_adf)
        frequencies = job_freq.frequencies
        
        print(run(frequencies))

Once you run the script an input file for the *ADF* and *Orca* jobs are
created. The *ADF* input looks like

::

    Atoms
          1         C      2.419290      0.606560      0.000000 
          2         C      1.671470      1.829570      0.000000 
          3         N      1.065290      2.809960      0.000000 
          4         H      2.000000      0.000000      1.000000 
          5         H      2.000000      0.000000     -1.000000 
          6         H      3.600000      0.800000      0.000000 
    End

    Basis
      Type DZP
    End

    Constraints
      Atom 1
      Atom 2
      Atom 3
    End

    Geometry
      Optim cartesian
    End

    Integration
      Accint 6.0
    End

    Scf
      Converge 1e-06
      Iterations 100
    End

 Running in a supercomputer 
----------------------------

Running in **Cartesius** or **Bazis** through the *Slurm* resource
manager can be done using and script like

.. code:: bash

    #!/bin/bash
    #SBATCH -t 00:10:00
    #SBATCH -N 1
    #SBATCH -n 8

    module load orca
    module load adf/2016.102

    source activate qmflows
    python optimization_ADF_freq_ORCA.py

The Slurm output looks like:

\`\`\` load orca/3.0.3 (PATH) discarding
/home/user/anaconda3/envs/qmflows/bin from PATH prepending
/home/user/anaconda3/envs/qmflows/bin to PATH [11:17:59] PLAMS working
folder: /nfs/home/user/orca/Opt/example/plams.23412 +-(running jobs) \|
Running adf ... [11:17:59] Job ADFjob started [11:18:18] Job ADFjob
finished with status 'successful' [11:18:18] Job ORCAjob started
[11:18:26] Job ORCAjob finished with status 'successful'

[ 0. 0. 0. 0. 0. 0. -360.547382 -360.14986 953.943089 954.3062 1049.2305
1385.756519 1399.961717 1399.979552 2602.599662 3080.45671 3175.710785
3177.612274] \`\`\`

 7. How the run function works? 
--------------------------------

 A little discussion about graphs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*qmflows* is meant to be used for both workflow generation and
execution. When you write a python script representing a workflow you
are explicitly declaring set of computations and their dependencies. For
instance the following workflow represent *ADF* and *Orca* computations
of the aforementioned example. In this
`graph <https://en.wikipedia.org/wiki/Graph_theory>`__ the octagons
represent quantum simulation using a package, while the ovals represent
both user input or data extracted from a simulation. Finally, the arrows
(called edges) represent the dependencies between all these objects.

**QMflows** automatically identify the dependencies between computations
and run them in the correct order (if possible in parallel).

 Restarting a simulation 
~~~~~~~~~~~~~~~~~~~~~~~~~

If you are running many computationally expensive calculations in a
supercomputer, it can happen that the computations take more time than
that allowed by the resource manager in your supercomputer and the
workflows gets cancel. But do not worry, you do not need to re-run all
the computations. Fortunately, *QMflows* offers a mechanism to restart
the workflow computations.

When running a workflow you will see that *QMflows* creates a set of
files called ``cache``. These files contain the information about the
workflow and its calculation. **In order to restart a workflow you only
need to relaunch it**, that's it!

 8. Advanced Examples 
----------------------

 Conditional Workflows 
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from noodles import gather
    from qmflows import dftb, adf, orca, run, Settings, templates, molkit, find_first_job

    # This examples illustrates the possibility to use different packages interchangeably.
    # Analytical frequencies are not available for B3LYP in ADF
    # This workflow captures the resulting error and submits the same job to ORCA.

    # Define the condition for a successful calculation
    def is_successful(result):
        return result.status not in ["failed", "crashed"]

    # Generate water molecule
    water = molkit.from_smiles('[OH2]', forcefield='mmff')

    # Pre-optimize the water molecule
    opt_water = dftb(
         templates.geometry, water, job_name="dftb_geometry")

    jobs = []

    # Generate freq jobs for 3 functionals
    for functional in ['pbe', 'b3lyp', 'blyp']:
        s=Settings()
        s.basis = 'DZ'
        s.functional = functional
        # Try to perform the jobs with adf or orca
        # take result from  first successful calculation
        freqjob = find_first_job(
              is_successful, [adf, orca], templates.freq.overlay(s), 
              opt_water.molecule, job_name=functional)
        jobs.append(freqjob)

    # Run workflow
    results = run(gather(*jobs), n_processes=1)

After running the above script you have a table like

::

    pbe     1533.267   3676.165   3817.097
    b3lyp   1515.799   3670.390   3825.813
    blyp    1529.691   3655.573   3794.110

 Non-adiabatic couplings 
~~~~~~~~~~~~~~~~~~~~~~~~~

`qmflows-namd <https://github.com/SCM-NV/qmflows-namd>`__ is a package
based on **QMflows** to compute the Non-adiabatic couplings for large
system involving thr use of **QMflows**, `Cython <http://cython.org/>`__
and `Numpy <http://www.numpy.org/>`__.

 9. Exception Handling 
-----------------------

Suppose you have a set of non-dependent calculations, for example single
point calculations coming from a molecular dynamic trajectory, as shown
in the figure below

If one of the single point calculations fails, the rest of the point in
the workflow will keep on running and the failed job will return a
**None** value for the requested property.

If the single point calculation would be the dependency of another
quantum calculation then the computation will crash.
