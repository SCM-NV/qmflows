Running a workflow
~~~~~~~~~~~~~~~~~~
A workflow in **Qmflows** consist of a set of computations and the dependencies between them,
explicitly declared by the user in the python script. This dependecies and the relation between
them form an graph (specifically an acyclic direct graph) that represent these relations.

**Qmflows** Builds the aforemention graph in order to realize the workflow evaluation order. For instance the figure below represent a simulation where firstly a molecular geometry optimization is carried out using the *ADF* package and some user defined ``Settings`` for the *ADF* simulation package. Subsequently, using the optimized molecular geometry from the previous step and another ``Settings`` for an orca simulation a job to compute the molecular frequencies is carried out.  

.. image:: _images/simple_graph.png

A python script corresponding with this graph can be::

   from plams import Molecule
   from qmflows import (adf, orca, run, Settings)

   # ADF optimization
   optmized_mol_adf = adf(inp, acetonitrile, job_name='acetonitrile_opt')

   # Orca Settings definition
   s2 = Settings()
   s2.specific.orca.main = "freq"
   s2.specific.orca.basis.basis = 'sto_sz'
   s2.specific.orca.method.functional = 'lda'
   s2.specific.orca.method.method = 'dft'

   # Orca Frequencies 
   job_freq = orca(s2, optmized_mol_adf)

   # Extract the frequencies from the Orca job
   frequencies = job_freq.frequencies

   # Run the graph
   result = run(frequencies)
   print(result)

   
Up to the invocation of the ``run`` function none of the computations have been executed,
it is the ``run`` function which builds and executes the dependencies. Since *Qmflows* needs to figure out
all the dependecies in the script, the ``run`` function takes as argument last dependency (or inner most dependy),
which in this case are the frequencies. The reason behind this, is that from the last dependency it is possible to
retrace all the dependecies.

**QWorks** uses  the noodles_ library under the hook to takes care of the construction and
execution of the dependecy graph.

.. _noodles: http://nlesc.github.io/noodles/

