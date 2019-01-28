General Overview & Getting Started
==================================

A basic recipe for running the quantum dot builder:

1.  Create two directories named ‘core’ and ‘ligand’. The 'core' directory should contain the input cores & the 'ligand' should contain the input ligands. The quantum dots will be exported to the 'QD' directory.

2. 	In qd_input.py (see qd-example_) alter the path variable. Path should point to the directory containing the directories mentioned in step 1.

3.  Enter all input cores and ligands as a list in the input_cores and input_ligands variables.

4.	Alter any optional argument in input_cores, input_ligands and/or argument_dict (see below).

5.	Start the job by running qd_input.py: \ ``python qd_input.py``

.. _qd-example: https://github.com/SCM-NV/qmflows/blob/master/test/QD_input_examples
