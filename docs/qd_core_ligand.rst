Input cores & ligands
=====================

Thia section related relates the importing and processing of cores and ligands. Ligand & cores can be imported from a wide range of different files and files types, which can roughly be divided into three categories:

1.  Files containing coordinates of a single molecule: .xyz, .pdb & .mol
2.  Python objects: ``plams.Molecule``, ``rdkit.Chem.Mol`` & (SMILES) ``str``
3.  Containers containing on or multiple input molecules: .txt, .xlsx & directories

In the later case, the container can consist of multiple SMILES strings or paths to .xyz, .pdb and/or .mol files. If necessary, containers are searched recursively. Both absolute and relative paths are explored.

Optional arguments
~~~~~~~~~~~~~~~~~~~

guess_bonds = False
-------------------

``Bool`` :
Try to guess bonds and bond orders in a molecule based on the types atoms and the relative of atoms.
Is set to False by default, with the exception of .xyz files.

column = 0
----------

``int`` :
The column containing the to be imported molecules.
Relevant when importing structures from .txt and .xlsx files with multiple columns.
Numbering starts from 0.

row = 0
-------

``int`` :
The first row in a column which contains a molecule.
Useful for when, for example, the very first row contains the title of aforementioned row, in which case row = 1 would be a sensible choice.
Relevant for .txt and .xlsx files.
Numbering starts from 0.

sheet_name = 'Sheet1'
---------------------

``str`` :
The name of the sheet containing the to be imported molecules.
Relevant for .xlsx files with multiple sheets or sheets with custom names.

core_indices = []
-----------------

``list`` consisting of ``int`` :
Manually specify the atomic index of one ore more atom(s) in the core that will be replaced with ligands. 
If left empty, all atoms of a user-specified element (see argument_dict: dummy) will be replaced with ligands.
The numbering atoms starts from 1, following the PLAMS 
[`1 <https://github.com/SCM-NV/PLAMS>`_,
`2 <https://www.scm.com/doc/plams/index.html>`_]
convention.

ligand_indices = []
-------------------

``tuple`` consisting of 1 or 2 ``int`` :
Manually specify the atomic index of the ligand atom that will be attached to core (implying argument_dict: split = False). 
If two atomic indices are provided, the bond between tuple[0] and tuple[1] will be broken and the molecule containing tuple[0] is attached to the core, (implying argument_dict: split = True).
Serves as an alternative to the functional group based substructure_search(), which identifies the to be atatched atom based on specific connectivity patterns (*i.e.* functional groups).
The numbering atoms starts from 1, following the PLAMS 
[`1 <https://github.com/SCM-NV/PLAMS>`_,
`2 <https://www.scm.com/doc/plams/index.html>`_]
convention.
