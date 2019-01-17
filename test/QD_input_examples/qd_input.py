import yaml
import qmflows.qd as QD
import os


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = '/Users/basvanbeek/Documents/CdSe/Week_5'

# The input cores from path/core/
input_cores = yaml.load("""
-   - Cd16Se13.xyz
    - guess_bonds: False
""")

# The input ligands from path/ligand/
input_ligands = yaml.load("""
- OCCCC
""")


# Optional arguments: these can safely be left to their default values
argument_dict = yaml.load("""
dir_name_list: [core, ligand, QD]
dummy: Cl
database_name: [ligand_database, QD_database]
use_database: True
ligand_opt: True
ligand_crs: True
qd_opt: False
qd_int: False
maxiter: 500
split: True
""")


# Runs the script: the ligand, core and quantum dot lists are returned
qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
ligand = ligand_list[0]
