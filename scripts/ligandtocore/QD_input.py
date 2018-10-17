import yaml
import QD


# Mandatory arguments: input_cores, input ligands & parhwill have to be specified by the user

# The location of the (to be created) working folder
path = r'/Users/basvanbeek/Documents/CdSe/Week_5'

# The input cores
input_cores = yaml.load("""
-   - Cd68Se55.xyz
    - guess_bonds: False
""")

# The input ligands
input_ligands = yaml.load("""
- OC
- OCC
- OCCC
- OCCCC
- OCCCCC
""")

"""
Input_cores & input_ligands:
A list of input cores and ligands.
Supported file types: .xyz, .pdb, .mol, folders, .txt, .xlsx, PLAMS molecules, RDKit molecules.
    .xyz:       Import an .xyz file. Guesses connectivity by default.
    .pdb:       Import a .pdb file.
    .mol:       Import a .mol file.
    folder:     Import all supported files from path/ligand/folder/. If None or '' is provided as.
                argument, all supported files from path/ligand/ will be imported.
    .txt:       Import all supported files from a .txt file (i.e. SMILES strings or
                                                             paths to .xyz files).
    .xlsx:      Import all supported files from an .xlsx file (i.e. SMILES strings or
                                                               paths to .xyz files).
    PLAMS mol:  Import a PLAMS Molecule object.
    RDKit mol:  Import a RDKit Chem.rdchem.Mol object.

Optional arguments:
    guess_bonds: bool
        Try to guess bonds in the molecule based on types and positions of atoms.
        By default this option is only enabled for .xyz files, as it is assumed that e.g. .pdb files
        already contain information about connectivity.

    column: int
        The column containing the to be imported molecules.
        Relevant for .txt and .xlsx files.

    row: int
        Ignore the first i rows within a column of to be imported molecules.
        Relevant for .txt and .xlsx files.

    sheet_name: str
        The name of the sheet containing the to be imported molecules
        Relevant for .xlsx files.
"""


# Optional arguments: these can safely be left to their default values
argument_dict = yaml.load("""
dir_name_list: [core, ligand, QD]
dummy: Cl
core_indices: []
ligand_indices: []
database_name: ligand_database.xlsx
use_database: True
core_opt: False
ligand_opt: True
qd_opt: False
maxiter: 10000
split: True
""")

"""
Optional arguments:
    dir_name_list <list>: [core, ligand, QD]
        A list containing the names of the be created directories.

    dummy <int> or <str>: Cl
        The atomic number (int) of atomic symbol (str) of the atoms in the core that should be
        replaced with ligands. Alternatively, dummy atoms can be manually specified with
        the core_indices argument.

    core_indices <list>[<int>]: []
        Specify the indices (int) of the core dummy atoms that should be replaced with ligands.
        Alternatively, all atoms of a given element can be identified as dummy atoms with the
        dummy argument.

    ligand_indices <list>[<int>] or <list>[<list>[<int>, <int>]] = []
        Specify the indices (int) of the ligand atoms that should be attached to the core.
        Alternatively, all atoms of a given element can be identified as dummy atoms with the
        dummy argument.

    database_name <str>: ligand_database.xlsx
        Name plus extension (str) of the (to be) created database where all results will be stored.
        When possible, ligands will be pulled from this database instead of
            reoptimizing their geometry.

    use_database <bool>: True
        Enables or disables the use of database_name.

    ligand_opt <False>: True
        Approach the global minimum (RDKit UFF) by systematically scanning dihedral angles.
        Requires well-defined bonds.

    core_opt <bool>: False
        Approach the global minimum (RDKit UFF) by systematically scanning dihedral angles.
        Requires well-defined bonds.

    qd_opt <bool>: False
        Optimize the quantum dot (qd, i.e core + all ligands) using ADF UFF.

    maxiter <int>:
        The maximum number of geometry iterations (int) used in qd_opt.

    split <bool>: True
        False: The ligand is the be attached to the core in its entirety.
            Examples:
                HO2CR -> HO2CR
                NH4+ -> NH4+
                -O2CR -> -O2CR
        True: A proton, counterion or functional group is to be removed from the ligand first.
            Examples:
                HO2CR -> -O2CR
                X-.NH4+ -> NH4+
                Na+.-O2CR -> -O2CR
                PhO2CR -> -O2CR
"""


# Runs the script: add ligand to core and optimize (UFF) the resulting qd with the core frozen
qd_list = QD.prep(input_ligands, input_cores, path, argument_dict)
