__all__ = ['prep']

import itertools
import time

from scm.plams import (Atom, MoleculeError, Settings)
import scm.plams.interfaces.molecule.rdkit as molkit

from .components import qd_functions as QD_scripts
from .components import qd_database as QD_database
from .components import qd_import_export as QD_inout
from .components import qd_ams as QD_ams


def prep(input_ligands, input_cores, path, arg):
    """
    function that handles all tasks related to prep_core, prep_ligand and prep_qd.

    input_ligands <list>[<plams.Molecule>]: A list of all input ligands.
    input_cores <list>[<plams.Molecule>]: A list of all input cores.
    path <str>: The path where all results will be stored.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A list of all quantum dots (core + n*ligands).
    """
    # The start
    time_start = time.time()
    print('\n')

    # Create the result directories (if they do not exist) and ligand and core lists
    folder_list = [QD_inout.create_dir(name, path) for name in arg['dir_name_list']]
    ligand_list = QD_inout.read_mol(input_ligands, folder_list[1])
    core_list = QD_inout.read_mol(input_cores, folder_list[0], is_core=True)

    # Adds the indices of the core dummy atoms to core.properties.core
    for core in core_list:
        prep_core(core, arg)

    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        database = QD_database.read_database(path, arg['database_name'])
    else:
        database = False

    # Optimize all ligands and find their functional groups
    ligand_list = list(prep_ligand(ligand, database, arg) for ligand in ligand_list)
    ligand_list = list(itertools.chain(*ligand_list))
    if not ligand_list:
        raise IndexError('No valid ligand functional groups were found, aborting run')

    # Write new entries to the ligand database
    if arg['use_database']:
        QD_database.write_database(ligand_list, database, path)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_list = list(prep_qd(core, ligand, folder_list[2]) for core in core_list for
                   ligand in ligand_list)

    # Optimize the qd with the core frozen
    if arg['qd_opt']:
        QD_ams.check_sys_var()
        if not qd_list:
            raise IndexError('No valid quantum dots were found, aborting geometry optimization')
        else:
            qd_list = list(QD_ams.ams_job(qd, maxiter=arg['maxiter'], job='qd_opt') for
                           qd in qd_list)

    # Calculate the (mean) interaction between ligands on the quantum dot surface
    if arg['qd_int']:
        if not qd_list:
            raise IndexError('No valid quantum dots were found, aborting single points')
        else:
            qd_list = list(QD_scripts.qd_int(qd) for qd in qd_list)
            if arg['use_database']:
                QD_database.write_database_qd(qd_list, path)

    # The End
    time_end = time.time()
    print('\nTotal elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')

    return qd_list, core_list, ligand_list


def prep_core(core, arg):
    """
    Function that handles all core operations.

    core <plams.Molecule>: The core molecule.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    dummy = arg['dummy']
    opt = arg['core_opt']

    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    if isinstance(dummy, str):
        dummy = Atom(symbol=dummy).atnum

    # Optimize the core with RDKit UFF if opt = True. Returns a RDKit molecule
    if opt:
        core = molkit.global_minimum(core)

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    if not core.properties.dummies:
        core.properties.dummies = [atom for atom in core.atoms if atom.atnum == dummy]
    else:
        core.properties.dummies = [core[index] for index in core.properties.dummies]
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Returns an error if no dummy atoms were found
    if not core.properties.dummies:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')


def prep_ligand(ligand, database, arg):
    """
    Function that handles all ligand operations,

    ligand <plams.Molecule>: The ligand molecule.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A copy of the ligand for each identified functional group.
    """
    opt = arg['ligand_opt']
    split = arg['split']
    crs = arg['ligand_crs']

    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand = QD_scripts.optimize_ligand(ligand, database, opt)

    # Identify functional groups within the ligand and add a dummy atom to the center of mass.
    if not ligand.properties.dummies:
        ligand_list = QD_scripts.find_substructure(ligand, split)
    else:
        if isinstance(ligand.properties.dummies, int):
            ligand.properties.dummies += -1
            split = False
        elif isinstance(ligand.properties.dummies, (list, tuple)):
            ligand.properties.dummies = [i - 1 for i in ligand.properties.dummies]
            split = True
        ligand_list = [QD_scripts.find_substructure_split(ligand, ligand.properties.dummies, split)]

    if crs:
        for ligand in ligand_list:
            QD_ams.ams_job(ligand, job='ligand_sp')

    return ligand_list


def prep_qd(core, ligand, qd_folder):
    """
    Function that handles all quantum dot (qd, i.e. core + all ligands) operations.

    core <plams.Molecule>: The core molecule.
    ligand <plams.Molecule>: The ligand molecule.
    qd_folder <str>: The quantum dot export folder.

    return <plams.Molecule>: The quantum dot (core + n*ligands).
    """
    # Rotate and translate all ligands to their position on the core.
    core_dummies = core.properties.dummies
    core = core.copy()
    core.properties.dummies = [core.closest_atom(dummy.coords) for dummy in core_dummies]
    ligand_list = [QD_scripts.rotate_ligand(core, ligand, i, core_dummy) for i, core_dummy in
                   enumerate(core.properties.dummies)]
    ligand_list, ligand_indices = zip(*ligand_list)
    core.delete_atom(core[-1])

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    qd = QD_scripts.combine_qd(core, ligand_list)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in ligand_indices]
    qd_indices += [i + 1 for i, atom in enumerate(core)]

    qd.properties = Settings()
    qd.properties.indices = qd_indices
    qd.properties.name = core.properties.name + '__' + ligand.properties.name
    qd.properties.path = qd_folder
    QD_inout.export_mol(qd, message='core + ligands:\t\t\t')

    return qd
