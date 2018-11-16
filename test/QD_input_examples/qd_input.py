import yaml
import qmflows.qd as QD
import os


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = os.get_cwd()

input_cores = yaml.load("""
-   - Cd16Se13.xyz
    - guess_bonds: False
""")

input_ligands = yaml.load("""
- OC
- OCC
- OCCC
- OCCCC
""")


# Optional arguments: these can safely be left to their default values
argument_dict = yaml.load("""
dir_name_list: [core, ligand, QD]
dummy: Cl
database_name: [ligand_database.xlsx, QD_database.xlsx]
use_database: True
core_opt: False
ligand_opt: True
ligand_crs: False
qd_opt: True
qd_int: True
maxiter: 900
split: True
""")


# Runs the script: the ligand, core and quantum dot lists are returned
qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
