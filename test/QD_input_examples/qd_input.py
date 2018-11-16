import yaml
import qmflows.qd as QD
import os


# Mandatory arguments: input_cores, input ligands & path will have to be specified by the user
path = r'/Users/basvanbeek/Documents/CdSe/Week_5'


input_cores = yaml.load("""
-   - Cd16Se13.xyz
    - guess_bonds: False
""")

# core_indices: [677, 757, 670, 771, 707, 759, 692, 765, 709, 687, 712, 671, 741, 723, 724, 683]

input_ligands = yaml.load("""
- input_ligands.txt
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


# Runs the script: add ligand to core and optimize (UFF) the resulting qd with the core frozen
qd_list, core_list, ligand_list = QD.prep(input_ligands, input_cores, path, argument_dict)
qd = qd_list[0]
ligand = ligand_list[0]
core = core_list[0]
