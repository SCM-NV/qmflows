=====
Step 4:		Optional arguments
=====

- dir_name_list = [core, ligand, QD]
:mod:`list` consisting of :mod:`str`: 
The names of the to be created folders.
dir_name_list[0] is the core input folder, dir_name_list[1] the ligand input folder and dir_name_list[2] the quantum dot output folder.

- dummy = Cl
:mod:`int`:
The atomic number  or atomic symbol :mod:`str` of the atoms in the core that should be replaced with ligands. Alternatively, dummy atoms can be manually specified with the core_indices variable.

- use_database = True
:mod:`bool`:
Enables or disables the use of database_name.

- ligand_opt = True
:mod:`bool`:
split the ligand into linear fragments and then recombine these fragments, searching for the optimal dihedral angle of the newly (re-)formed bond in the process. Involved an optimization with RDKit UFF.

- split = True
:mod:`bool`:
If False: The ligand is to be attached to the core in its entirety.

    NR\ :sub:`4`\ :sup:`+` \     ->     NR\ :sub:`4`\ :sup:`+` \
    
    O\ :sub:`2`\CR        ->     O\ :sub:`2`\CR
    
    HO\ :sub:`2`\CR       ->     HO\ :sub:`2`\CR
    
    H\ :sub:`3`\CO\ :sub:`2`\CR     ->     H\ :sub:`3`\CO\ :sub:`2`\CR

If True: A proton, counterion or functional group first is to be removed first from the ligand.

    X\ :sup:`-`\.NR\ :sub:`4`\     ->     NR\ :sub:`4`\ :sup:`+` \
    
    HO\ :sub:`2`\CR       ->     O\ :sup:`-`\ :sub:`2`\CR
    
    Na\ :sup:`+`\.O\ :sup:`-`\ :sub:`2`\CR	-> 	O\ :sup:`-`\ :sub:`2`\CR
    
    H\ :sub:`3`\CO\ :sub:`2`\CR     ->     O\ :sup:`-`\ :sub:`2`\CR

- ligand_crs = False
:mod:`bool`:
Calculate the ligand volume, surface area and octanol/water partition coefficient
with ADF MOPAC + COSMO-RS.

- qd_opt = False
:mod:`bool`:
Optimize the quantum dot (i.e. core + all ligands) with ADF UFF.
The geometry of the core and ligand atoms directly attached to the core are frozen during this optimization.

- maxiter = 500
:mod:`int`:
The maximum number of geometry iterations during qd_opt.

- qd_int = False
:mod:`bool`:
Perform an activation strain analyses on the ligands attached to the quantum dot surface with RDKit UFF. The core is removed during this process; the analyses is thus exclusively focused on ligand deformation and inter-ligand interaction.
Yields three terms:

d\ *E*\ :sub:`strain`\ : 	The energy required to deform the ligands from their equilibrium geometry to the geometry they adopt on the quantum dot surface. This term is, by definition, destabilizing.

d\ *E*\ :sub:`int`\ :	The mutual interaction between all deformed ligands. This term is characterized by the non-covalent interaction between ligands (UFF Lennard-Jones potential) and, depending on the inter-ligand distances, can be either stabilizing or destabilizing.

d\ *E* :	The sum of d\ *E*\ :sub:`strain`\  and d\ *E*\ :sub:`int`\  accounts for both the destabilizing ligand deformation and (de-)stabilizing interaction between all ligands in the absence of the core.
