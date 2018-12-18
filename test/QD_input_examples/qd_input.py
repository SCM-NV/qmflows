import yaml
import qmflows.qd as QD
import os


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = '/Users/basvanbeek/Documents/CdSe/Week_5'

# The input cores from path/core/
input_cores = yaml.load("""
-   - Cd176Se147.xyz
    - guess_bonds: False
    - core_indices: [338, 375, 327, 372, 341, 355, 374, 356, 340]
""")

# The input ligands from path/ligand/
input_ligands = yaml.load("""
- OC(CCCCCC)=O
""")


# Optional arguments: these can safely be left to their default values
argument_dict = yaml.load("""
dir_name_list: [core, ligand, QD]
dummy: Cl
database_name: [ligand_database, QD_database]
use_database: True
ligand_opt: True
ligand_crs: False
qd_opt: True
qd_int: False
qd_dissociate: True
maxiter: 500
split: True
""")


# Runs the script: the ligand, core and quantum dot lists are returned
qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
