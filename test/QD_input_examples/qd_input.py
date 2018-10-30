import yaml
import qmflows.qd as QD
import os


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = r'/Users/bvanbeek/Documents/CdSe/Week_5'

input_cores = yaml.load("""
-   - Cd68Se55.xyz
    - guess_bonds: False
""")

input_ligands = yaml.load("""
- OCCCCC
""")

"""
path <str>:
    The location of the (to be created) working folders <str>.
    Set to os.getcwd() to use the current directory.

<mol> = <str>, <scm.plams.core.basemol.Molecule> or <rdkit.Chem.rdchem.Mol>
Input_cores <list>[<mol>] or <list>[<list>[<mol>, <dict>]]:
Input_ligands <list>[<mol>] or <list>[<list>[<mol>, <dict>]]:
    A list of all input molecules <mol>.
    The script will automatically determine the supplied file type.
    Supported file types:
        .xyz:       Import an .xyz file. Guesses connectivity by default.
        .pdb:       Import a .pdb file.
        .mol:       Import a .mol file.
        folder:     Import all supported files from path/ligand/folder/.
                    If None or '' is provided as argument,
                    all supported files from path/ligand/ will be imported instead.
        .txt:       Import all supported files from a .txt file (e.g. SMILES strings or
                                                                 paths to .xyz files).
        .xlsx:      Import all supported files from an .xlsx file (e.g. SMILES strings or
                                                                   paths to .xyz files).
        PLAMS mol:  Import a PLAMS Molecule object.
        RDKit mol:  Import a RDKit Chem.rdchem.Mol object.

    To supply optional arguments <dict>, <mol> has to be exchanged for <list>[<mol>, <dict>].
    Optional arguments:
        guess_bonds <bool>: False
            Try to guess bonds in the molecule based on types and positions of atoms.
            Is set to False by default, with the exception of .xyz files.

        column <int>: 0
            The column containing the to be imported molecules.
            Relevant for .txt and .xlsx files.

        row <int>: 0
            Ignore the first i rows within a column of to be imported molecules.
            Relevant for .txt and .xlsx files.

        sheet_name <str>: Sheet1
            The name of the sheet containing the to be imported molecules
            Relevant for .xlsx files.

        core_indices <list>[<int>]: []
            Manually specify the index <int> of the core atoms that should be replaced with ligands.
            Alternatively, all atoms of a given element can be identified as dummy atoms with the
            dummy argument (see argument_dict).

        ligand_indices <list>[<int>] or <list>[<int>, ..., <int>]: []
            Manually specify the index <int> of the ligand atom that should be attached to the core.
            Serves as an alternative to the functional group based substructure_search().
            If a list of ints is supplied instead, list[0] will be attached to the core and list[-1]
                will be removed from the ligand (see the keyword 'split' in argument_dict).

"""


# Optional arguments: these can safely be left to their default values
argument_dict = yaml.load("""
dir_name_list: [core, ligand, QD]
dummy: Cl
database_name: ligand_database.xlsx
use_database: True
core_opt: False
ligand_opt: True
ligand_crs: False
qd_opt: False
qd_int: True
maxiter: 2000
split: True
""")

"""
Optional arguments:
    dir_name_list <list>[<str>]: [core, ligand, QD]
        A list containing the names of the be created directories.

    dummy <int> or <str>: Cl
        The atomic number <int> or atomic symbol <str> of the atoms in the core that should be
        replaced with ligands. Alternatively, dummy atoms can be manually specified with
        the core_indices argument.

    database_name <str>: ligand_database.xlsx
        Name plus extension <str> of the (to be) created database where all results will be stored.
        When possible, ligands will be pulled from this database instead of
            reoptimizing their geometry.

    use_database <bool>: True
        Enables or disables the use of database_name.

    ligand_opt <bool>: True
    core_opt <bool>: False
        Approach the global minimum (RDKit UFF) by systematically scanning dihedral angles.
        Requires well-defined bonds.
        WARNING: No well-defined bonds are present in the cores be default, enable at your own risk!

    ligand_crs <bool>: False
        Calculate various ligand properties using ADF MOPAC + COSMO-RS.

    qd_opt <bool>: False
        Optimize the quantum dot (qd, i.e core + all ligands) using ADF UFF.

    qd_int <bool>: False
        Perform an activation strain analyses on the ligands attached to the quantum dot surface
            using RDKit UFF.

    maxiter <int>: 2000
        Maxiter + 50 is the maximum number of geometry iterations <int> used in qd_opt.

    split <bool>: True
        False: The ligand is to be attached to the core in its entirety.
            Examples:
                NR4+      -> NR4+
                -O2CR     -> -O2CR
                HO2CR     -> HO2CR
                PhO2CR    -> PhO2CR
        True: A proton, counterion or functional group first has to be removed from the ligand.
            Examples:
                X-.NR4+   -> NR4+
                HO2CR     -> -O2CR
                Na+.-O2CR -> -O2CR
                PhO2CR    -> -O2CR
"""


# Runs the script: add ligand to core and optimize (UFF) the resulting qd with the core frozen
qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
qd = qd_list[0]
ligand = ligand_list[0]
core = core_list[0]
