import copy
import itertools
import time

from scm.plams import (Atom, MoleculeError, Settings)
import scm.plams.interfaces.molecule.rdkit as molkit

import qd_functions as QD_scripts
import qd_database as QD_database
import qd_import_export as QD_inout


def prep(input_ligands, input_cores, path, arg):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.
    """
    # The start
    time_start = time.time()
    print('\n')

    # Create the result directories (if they do not exist), ligands and cores
    folder_list = [QD_inout.create_dir(name, path) for name in arg['dir_name_list']]
    ligand_list = QD_inout.read_mol(input_ligands, folder_list[1])
    core_list = QD_inout.read_mol(input_cores, folder_list[0], is_core=True)

    # Adds the indices of the core dummy atoms to core.properties.core
    for core in core_list:
        prep_core(core, arg['core_indices'], arg['dummy'], arg['core_opt'])

    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        database = QD_database.read(folder_list[1], arg['database_name'])
    else:
        database = False

    # Optimize all ligands and find their functional groups
    ligand_list = list(prep_ligand(ligand, database, arg['ligand_indices'], arg['ligand_opt'],
                                   arg['split']) for ligand in ligand_list)
    ligand_list = list(itertools.chain(*ligand_list))
    if not ligand_list:
        raise IndexError('No valid ligand functional groups were found, aborting run')

    # Write new entries to the ligand database
    if arg['use_database']:
        QD_database.write(ligand_list, database)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_list = list(prep_qd(core, ligand, folder_list[2]) for core in core_list for
                   ligand in ligand_list)

    # Check if the ADF environment variables are set and optimize the qd with the core frozen
    if arg['qd_opt']:
        QD_scripts.check_sys_var()
        if not qd_list:
            raise IndexError('No valid quantum dots were found, aborting geometry optimization')
        for qd in qd_list:
            QD_scripts.ams_job(qd, arg['maxiter'])

    # The End
    time_end = time.time()
    print('\nTotal elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')

    return qd_list


def prep_core(core, core_dummies, dummy=0, opt=False):
    """
    Function that handles all core operations.
    """
    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum

    # Optimize the core with RDKit UFF if opt = True. Returns a RDKit molecule
    if opt:
        core = molkit.global_minimum(core)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    if not core_dummies:
        core_dummies = [atom for atom in core.atoms if atom.atnum == dummy]
    else:
        core_dummies = [core[index] for index in core_dummies]
    core.properties.core_dummies = core_dummies
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Returns an error if no dummy atoms were found
    if not core_dummies:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')


def prep_ligand(ligand, database, ligand_indices, opt=True, split=True):
    """
    Function that handles all ligand operations,
    """
    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand = QD_scripts.optimize_ligand(ligand, opt, database)

    # Identify functional groups within the ligand and add a dummy atom to the center of mass.
    if not ligand_indices:
        ligand_list = QD_scripts.find_substructure(ligand, split)
    else:
        if isinstance(ligand_indices, int):
            ligand_indices += -1
            split = False
        elif isinstance(ligand_indices, list):
            ligand_indices = [i - 1 for i in ligand_indices]
            split = True
        ligand_list = [QD_scripts.find_substructure_split(ligand, ligand_indices, split)]

    return ligand_list


def prep_qd(core, ligand, qd_folder):
    """
    Function that handles all quantum dot (qd, i.e. core + all ligands) operations.
    """
    # Rotate and translate all ligands to their position on the core.
    core = copy.deepcopy(core)
    ligand_list = [QD_scripts.rotate_ligand(core, ligand, i, core_dummy) for i, core_dummy in
                   enumerate(core.properties.core_dummies)]
    ligand_list, ligand_indices = zip(*ligand_list)
    core.delete_atom(core[-1])

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    qd = QD_scripts.combine_qd(core, ligand_list)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in ligand_indices]
    qd_indices += [i + 1 for i, atom in enumerate(core)]

    qd.properties = Settings()
    qd.properties.qd_indices = qd_indices
    qd.properties.name = core.properties.name + '__' + ligand.properties.name
    qd.properties.source_folder = qd_folder
    QD_inout.export_mol(qd, message='core + ligands:\t\t\t')

    return qd
