__all__ = ['prep']

import itertools
import time
import os
import pandas as pd

from scm.plams import (Atom, MoleculeError, Settings)

from .components import qd_functions as QD_scripts
from .components import qd_database as QD_database
from .components import qd_import_export as QD_inout
from .components import qd_ams as QD_ams
from .components import qd_dissociate as QD_dissociate
from .components import qd_ligand_opt as QD_ligand_opt


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

    # Raises an error if mol_list is empty
    if not ligand_list:
        raise IndexError('No valid input ligands were found, aborting run')
    elif not core_list:
        raise IndexError('No valid input cores were found, aborting run')

    # Adds the indices of the core dummy atoms to core.properties.core
    for core in core_list:
        prep_core(core, arg)

    # Optimize the ligands, find functional groups, calculate properties and read/write the results
    ligand_list = prep_ligand_1(ligand_list, path, arg)

    # Combine the core with the ligands, yielding qd, and format the resulting list
    qd_list = list(prep_qd_1(core, ligand, folder_list[2]) for core in core_list for
                   ligand in ligand_list)

    # Optimize the quantum dots, perform an activation strain analyses and read/write the results
    qd_list = prep_qd_2(qd_list, path, arg)

    # The End
    time_end = time.time()
    print('\n' + QD_scripts.get_time(),
          'Total elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec')

    return qd_list, core_list, ligand_list


def prep_core(core, arg):
    """
    Function that handles the identification and marking of all core dummy atoms.

    core <plams.Molecule>: The core molecule.
    arg <dict>: A dictionary containing all (optional) arguments.
    """
    # Checks the if the dummy is a string (atomic symbol) or integer (atomic number)
    dummy = QD_scripts.to_atnum(arg['dummy'])

    # Returns the indices (integer) of all dummy atom ligand placeholders in the core
    # An additional dummy atom is added at the core center of mass for orientating the ligands
    if not core.properties.dummies:
        core.properties.dummies = [atom for atom in reversed(core.atoms) if atom.atnum == dummy]
    else:
        core.properties.dummies.sort(reverse=True)
        core.properties.dummies = [core[index] for index in core.properties.dummies]
    core.add_atom(Atom(atnum=0, coords=(core.get_center_of_mass())))

    # Returns an error if no dummy atoms were found
    if not core.properties.dummies:
        raise MoleculeError(Atom(atnum=dummy).symbol +
                            ' was specified as dummy atom, yet no dummy atoms were found')


def prep_ligand_1(ligand_list, path, arg):
    """
    Function that handles ligand operations.
    Read/write the results from/to a ligand database, launch prep_ligand_2() and
        calculate properties with MOPAC + COSMO-RS.

    ligand_list <list>[<plams.Molecule>]: A list of all ligand molecules.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A copy of all ligands for each identified functional group.
    """
    # Open the ligand database and check if the specified ligand(s) is already present
    if arg['use_database']:
        ligand_database = QD_database.read_database(path, database_name='Ligand_database')
    else:
        ligand_database = None

    # Optimize all ligands and find their functional groups
    ligand_list = list(itertools.chain.from_iterable(prep_ligand_2(ligand, ligand_database, arg) for
                                                     ligand in ligand_list))
    if not ligand_list:
        raise IndexError('No valid ligand functional groups found, aborting run')

    if arg['ligand_crs']:
        QD_ams.check_sys_var()
        for ligand in ligand_list:
            QD_ams.ams_job_mopac_sp(ligand)

    # Write new entries to the ligand database
    if arg['use_database']:
        if not arg['ligand_opt']:
            for ligand in ligand_list:
                ligand.properties.entry = True
        QD_database.write_database(ligand_list, ligand_database, path, mol_type='ligand')

    return ligand_list


def prep_ligand_2(ligand, database, arg):
    """
    Function that handles ligand operations.
    Optimize the ligand and search for viable user-defined functional groups.

    ligand <plams.Molecule>: The ligand molecule.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A copy of the ligand for each identified functional group.
    """
    split = arg['split']

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

    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand_list = [QD_ligand_opt.optimize_ligand(ligand, database, arg['ligand_opt']) for
                   ligand in ligand_list if ligand_list]

    return ligand_list


def prep_qd_1(core, ligand, qd_folder):
    """
    Function that handles quantum dot (qd, i.e. core + all ligands) operations.
    Combine the core and ligands and assign properties to the quantom dot.

    core <plams.Molecule>: The core molecule.
    ligand <plams.Molecule>: The ligand molecule.
    qd_folder <str>: The quantum dot export folder.

    return <plams.Molecule>: The quantum dot (core + n*ligands).
    """
    # Rotate and translate all ligands to their position on the core.
    indices = [(dummy.get_atom_index(), -1, ligand.properties.dummies.get_atom_index(), -1) for
               i, dummy in enumerate(core.properties.dummies)]
    core = core.copy()
    ligand.add_atom(Atom(atnum=0, coords=ligand.get_center_of_mass()))
    ligand_list = [QD_scripts.rotate_ligand(core, ligand.copy(), atoms, residue_number=i+1) for
                   i, atoms in enumerate(indices)]
    ligand_atoms = [ligand.properties.anchor for ligand in ligand_list]

    # Attach the rotated ligands to the core, returning the resulting strucutre (PLAMS Molecule).
    core.merge_mol(ligand_list)
    qd = core
    for atoms in indices:
        qd.delete_atom(qd[atoms[0]])
    for atom in reversed(qd.atoms):
        if atom.atnum == 0:
            qd.delete_atom(atom)

    # indices of all the atoms in the core and the ligand heteroatom anchor.
    qd_indices = [qd.atoms.index(atom) + 1 for atom in qd if
                  atom.properties.pdb_info.ResidueName == 'COR' or atom in ligand_atoms]

    qd_name = core.properties.name + '__'
    qd_name += str(len(ligand_list)) + '_' + ligand.properties.name + '@' + ligand.properties.group

    qd.properties = Settings()
    qd.properties.indices = qd_indices
    qd.properties.name = qd_name
    qd.properties.path = qd_folder
    QD_inout.export_mol(qd, message='core + ligands:\t\t')

    return qd


def prep_qd_2(qd_list, path, arg):
    """
    Function that handles quantum dot (qd, i.e. core + all ligands) operations.
    Optimize the quantum dot, perform and activation strain analyses on the ligands and read/write
        the results from/to a quantom dot database.

    ligand_list <list>[<plams.Molecule>]: A list of all quantom dots.
    database <pd.DataFrame>: Database of previous calculations.
    arg <dict>: A dictionary containing all (optional) arguments.

    return <list>[<plams.Molecule>]: A list of all optimized quantom dots.
    """
    if not qd_list:
        raise IndexError('No valid quantum dots found, aborting')

    # Open the quantum dot database and check if the specified quantum dot(s) is already present
    if arg['use_database']:
        qd_database = QD_database.read_database(path, database_name='QD_database')
    else:
        qd_database = None

    # Optimize the qd with the core frozen
    if arg['qd_opt']:
        QD_ams.check_sys_var()
        qd_list = list(QD_ams.qd_opt(qd, qd_database, arg) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface
    if arg['qd_int']:
        print(QD_scripts.get_time(), 'calculating ligand distortion and inter-ligand interaction...')
        qd_list = list(QD_scripts.qd_int(qd) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of
    # one or more ligands
    if arg['qd_dissociate']:
        print(QD_scripts.get_time(), 'calculating ligand dissociation energy...')

        def diss_list_to_pd(diss_list, residue_list, top_dict):
            gen = ((tuple(res), tuple(top_dict[i] for i in res),
                   qd.properties.Eint, qd.properties.Estrain, qd.properties.E) for
                   qd, res in zip(diss_list, residue_list))
            keys = ('Residue numbers', 'Topology', 'Eint', 'Estrain', 'E')
            return pd.DataFrame(dict(zip(keys, zip(*gen))))

        for qd in qd_list:
            top_dict = QD_dissociate.get_topology_dict(qd, dist=4.5)
            qd_gen = QD_dissociate.dissociate_ligand_cd(qd)
            # entries = QD_dissociate.diss_list_to_pd(diss_list, residue_list, top_dict)
            # entries.to_excel(os.path.join(path, 'dissociate.xlsx'))

    # Write the new quantum dot results to the quantum dot database
    if arg['use_database']:
        if not arg['qd_opt']:
            for qd in qd_list:
                qd.properties.entry = True
        QD_database.write_database(qd_list, qd_database, path, mol_type='qd')

    return qd_list
