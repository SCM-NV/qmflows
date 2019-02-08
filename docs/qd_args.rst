Further Optional Arguments
==========================

There are a number of arguments which can be used to modify the functionality and behaviour of the quantum dot builder. Herein an overview is provided.

Optional Arguments
~~~~~~~~~~~~~~~~~~

dir_name_list = ['core', 'ligand', 'QD']
----------------------------------------

``list`` consisting of ``str`` : 
The names of the (to be created) folders.
By default, ligand structures will be stored and read from dir_name_list[0], cores will be stored and read dir_name_list[1] and the combined core+ligands will be stored and read from dir_name_list[2].
Structures can be read from different folders if their filename is prepended with its absolute path, thus effectively ignoring this argument.

dummy = 'Cl'
------------

``int`` or ``str`` :
The atomic number or atomic symbol of the atoms in the core that is to be replaced with ligands. 
Alternatively, dummy atoms can be manually specified with the core_indices variable.

use_database = True
-------------------

``bool`` :
Enables or disables the storing and pulling of structures and properties from a user-created database (stored in .json and .xlsx formats). The script will attempt to pull a structure from the database if a match is found between a current input ligand and/or core+ligands and a previously optimized structure.

ligand_opt = True
-----------------

``bool`` :
Optimize the geometry of the to be attached ligands. 
The ligand is split into one or multiple (more or less) linear fragments, which are subsequently optimized (RDKit UFF 
[`1 <http://www.rdkit.org>`_,
`2 <https://github.com/rdkit/rdkit>`_,
`3 <https://doi.org/10.1021/ja00051a040>`_]
) and reassembled while checking for the optimal dihedral angle. The ligand fragments are biased towards more linear conformations to minimize inter-ligand repulsion once the ligands are attached to the core.

split = True
------------

``bool`` :
If False: The ligand in its entirety is to be attached to the core.

    NR\ :sub:`4`\ :sup:`+` \                    ->     NR\ :sub:`4`\ :sup:`+` \
    
    O\ :sub:`2`\CR                              ->     O\ :sub:`2`\CR
    
    HO\ :sub:`2`\CR                             ->     HO\ :sub:`2`\CR
    
    H\ :sub:`3`\CO\ :sub:`2`\CR                 ->     H\ :sub:`3`\CO\ :sub:`2`\CR

If True: A proton, counterion or functional group is to be removed from the ligand before attachment to the core.

    X\ :sup:`-`\.NR\ :sub:`4`\                  ->     NR\ :sub:`4`\ :sup:`+` \
    
    HO\ :sub:`2`\CR                             ->     O\ :sup:`-`\ :sub:`2`\CR
    
    Na\ :sup:`+`\.O\ :sup:`-`\ :sub:`2`\CR	    -> 	   O\ :sup:`-`\ :sub:`2`\CR
    
    H\ :sub:`3`\CO\ :sub:`2`\CR                 ->     O\ :sup:`-`\ :sub:`2`\CR

ligand_crs = False
------------------

``bool`` :
Perform a property calculation with COSMO-RS 
[`4 <https://www.scm.com/doc/COSMO-RS/index.html>`_,
`5 <https://doi.org/10.1021/j100007a062>`_, 
`6 <https://doi.org/10.1021/jp980017s>`_, 
`7 <https://doi.org/10.1139/V09-008>`_]
; the COSMO surfaces are constructed using ADF MOPAC
[`8 <https://www.scm.com/doc/MOPAC/Introduction.html>`_, 
`9 <http://openmopac.net/>`_, 
`10 <https://doi.org/10.1007/s00894-012-1667-x>`_]
.
The following properties are calculated:
    
1. The surface area of the ligand (A\ :sup:`2`\) as defined by its COSMO surface.
    
2. The volume of the ligand (A\ :sup:`3`\) as defined by the volume encompassed by its COSMO surface.
    
3. The solvation energy of the ligand (kcal mol\ :sup:`-1`\), at infinite dilution, in the following solvents: acetone, acetonitrile, dimethyl formamide (DMF), dimethyl sulfoxide (DMSO), ethyl acetate, ethanol, *n*-hexane, toluene and water.

qd_opt = False
--------------

``bool`` :
Optimize the quantum dot (i.e. core + all ligands) with ADF UFF
[`3 <https://doi.org/10.1021/ja00051a040>`_,
`11 <https://www.scm.com/doc/UFF/index.html>`_]
.
The geometry of the core and ligand atoms directly attached to the core are frozen during this optimization.

maxiter = 500
-------------

``int`` :
The maximum number of iterations during the geometry optimization of the quantum dot.
Only applicable if qd_opt = True.

qd_int = False
--------------

``bool`` :
Perform an activation strain analyses
[`12 <https://doi.org/10.1002/9780470125922.ch1>`_,
`13 <https://doi.org/10.1002/wcms.1221>`_,
`14 <https://doi.org/10.1021/acs.jpcc.5b02987>`_] (kcal mol\ :sup:`-1`\)
on the ligands attached to the quantum dot surface with RDKit UFF
[`1 <http://www.rdkit.org>`_,
`2 <https://github.com/rdkit/rdkit>`_,
`3 <https://doi.org/10.1021/ja00051a040>`_]
. 
The core is removed during this process; the analyses is thus exclusively focused on ligand deformation and inter-ligand interaction.
Yields three terms:

1.  d\ *E*\ :sub:`strain`\  : 	The energy required to deform the ligands from their equilibrium geometry to the geometry they adopt on the quantum dot surface. This term is, by definition, destabilizing. Also known as the preperation energy (d\ *E*\ :sub:`prep`\).

2.  d\ *E*\ :sub:`int`\  :	The mutual interaction between all deformed ligands. This term is characterized by the non-covalent interaction between ligands (UFF Lennard-Jones potential) and, depending on the inter-ligand distances, can be either stabilizing or destabilizing.

3.  d\ *E* :	The sum of d\ *E*\ :sub:`strain`\  and d\ *E*\ :sub:`int`\ . Accounts for both the destabilizing ligand deformation and (de-)stabilizing interaction between all ligands in the absence of the core.
