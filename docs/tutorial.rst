What is QMWorks?
===============
Research on modern computational quantum chemistry relies on a set of computational tools to carry out calculations. The complexity of the calculations usually requires
intercommunication between the aforementioned tools, such communication is usually done through scripts  (e.g. Bash, Python, Perl, etc.) that try to automatize input/output actions like: launching the computations in a cluster, reading the resulting output and feeding the relevant numerical result to another program. Such scripts are difficult to maintain and extend, requiring a significant programming expertise to work with them. Being then desirable a set of automatic and extensible tools that allows to perfom complex simulations in heterogeneous hardware plaftorms. This library tackles the construction and efficient execution of computational chemistry workflows.

This allows computational chemists to use the emerging massively parallel compute environments in an easy manner and focus on interpretation of scientific data rather than on tedious job submission procedures and manual data processing. In this context, **QMWorks** stands as a robust and flexible python library, that allows both begginers and experts to automate the following tasks:

* Input generation (*Plams* and *Templates*)
* Dependencies management (*Noodles*).
* Effcient data storage (*HDF5*).
* Distribution in heterogeneous hardware platforms.
* Jobs failure detection and recovery.
* Control over the granularity of parallelism.

.. image:: _images/mindmap.png
    :scale: 75 %
  
  
Tutorial
========
In this tutorial we will walk through some minimal examples that show the **QWORKS** Capabilities and philosophy. like single point calculations, geometry optimizatios,
transition state optimizations, etc. Before we starts with computational chemistry calculation we first need to review the Templates concept.

Note:
  The first noteworthy features of *QMWorks* is that it uses *Python3.5* if you have previous experience with *Python2.* and want to know what are the mayor diferences
  have a look at python2andpython3_.


Working with templates
~~~~~~~~~~~~~~~~~~~~~~
**QMworks** Offers a :ref:`templates` set that contains some default values for the most frequent calculation in computational chemistry: single point, geometry optimization, TS optimization, frequency calculations, etc.In order to understand the logic behind the templates, suppose that you want to carry out a calculation, like an optimization follow by frequencies calculation.  Also, you want to use a predifined generic template that contains details that some defaults like numerical tolerance, number of scf cycles, etc.  In the figure below you can see an schematic representation of the user input and the generic template. *QMWorks* will take your input and overlay it with the template that you have pick from the default :ref:`templates`.

.. image:: _images/default_tree.jpg

This overlay will merge your keywords with the template keywords in such a way that if a property or method is present in either one of the two set the property or method will be part of the final input. On the other hand if a property or method is presented in both input sets  *the user keyword always has preference over the template* as depicted in the figure below
	   
.. image:: _images/merged_tree.jpg

The input keywords are specified using the |Settings| class. |Settings| is a :mod:`dict` subclass that represent the data in a hierarchical tree-like structure. for example ::

  from qmworks import Settings, templates

  s = Settings()  # (1)
  # generic keyword 
  s.basis = "DZP"  #  (2)
  # "specific" allows the user to apply specific keywords for a package
  s.specific.adf.basis.core = "large"  # (3)

  input = templates.singlepoint.overlay(s)  # (4)
  
The above code snippet shows how to create a |Settings| instance object (1), then the generic keyword *basis*  declares that the "DZP" should be used together with the *large* keyword
of **ADF**. Finally the user's keywords are merged with the defaults resultin in a input like ::

  Basis
    Core large
    Type DZP
  End

  Integration
    Accint 4.0
  End

  Scf
    Converge 1e-06
    Iterations 100
  End

  Xc
   Lda
  End

  
A simple Example
----------------
Suppose you have a molecule and you first want to calculate a first approximation to the structure of minimal energy
(e.g. using), then using this first approximation you want to perform a more accurate calculation for the optimized
structure.

	   
	   
Running
~~~~~~~


A more Advance Example
----------------------


.. _python2andpython3: https://wiki.python.org/moin/Python2orPython3

