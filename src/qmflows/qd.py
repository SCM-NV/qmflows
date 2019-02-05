__all__ = ['prep']

import itertools
import time
from os.path import join

from scm.plams import (Atom, MoleculeError)

from .components import qd_ams as QD_ams
from .components import qd_bde as QD_BDE
from .components import qd_database as QD_database
from .components import qd_functions as QD_scripts
from .components import qd_import_export as QD_inout
from .components import qd_ligand_opt as QD_ligand_opt
from .components import qd_ligand_rotate as QD_ligand_rotate


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
    cor_dir, lig_dir, qd_dir = [QD_inout.create_dir(name, path) for name in arg['dir_name_list']]
    ligand_list = QD_inout.read_mol(input_ligands, lig_dir)
    core_list = QD_inout.read_mol(input_cores, cor_dir, is_core=True)

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
    qd_list = list(QD_ligand_rotate.ligand_to_qd(core, ligand, qd_dir) for core
                   in core_list for ligand in ligand_list)

    # Optimize the quantum dots, perform an activation strain analyses and read/write the results
    qd_list = prep_qd(qd_list, path, arg)

    # The End
    time_end = time.time()
    message = '\n' + QD_scripts.get_time()
    message += 'Total elapsed time:\t\t' + '%.4f' % (time_end - time_start) + ' sec'
    print(message)

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
        core.properties.dummies = [atom for atom in core.atoms if atom.atnum == dummy]
    else:
         core.properties.dummies = [core[index] for index in core.properties.dummies]

    # Delete all core dummy atoms
    for at in reversed(core.properties.dummies):
        core.delete_atom(at)

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
            QD_ams.ams_job_mopac_crs(ligand)

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
        if len(ligand.properties.dummies) == 1:
            ligand.properties.dummies = ligand.properties.dummies[0] -1
            split = False
        elif len(ligand.properties.dummies) == 2:
            ligand.properties.dummies = [i - 1 for i in ligand.properties.dummies]
            split = True
        ligand_list = [QD_scripts.find_substructure_split(ligand, ligand.properties.dummies, split)]

    # Handles all interaction between the database, the ligand and the ligand optimization
    ligand_list = [QD_ligand_opt.optimize_ligand(ligand, database, arg['ligand_opt']) for
                   ligand in ligand_list if ligand_list]

    return ligand_list


def get_job_settings(arg_dict, jobs=1):
    """
    """
    if isinstance(arg_dict, bool):
        ret = [None for i in range(jobs*2)]
        print(QD_scripts.get_time() + 'No user-specified jobs & settings found for qd_dissociate, \
              switching to defaults')
    else:
        try:
            ret = [arg_dict[item] for item in arg_dict]
            len_ret = len(ret)
        except TypeError:
            raise TypeError('Only booleans, dictiories or dictionary derived objects are \
                            valid when defining jobs')

        # Pad with <None> ret if is smaller than 2 * *jobs*
        if len_ret < jobs*2:
            for i in range(jobs*2 - len_ret):
                ret.append(None)
            print(QD_scripts.get_time() + 'No jobs & settings have been specified found for the \
                  last ' + str(jobs - len_ret/2) + ' jobs, switching to defaults')
        # Pop entries from ret if it is larger than 2 * *jobs*
        elif len_ret > jobs*2:
            ret = ret[0:jobs*2]
            print(QD_scripts.get_time() + str(len_ret / 2) + ' jobs have been specified while the \
                  argument only support ' + str(jobs) + ', the last ' + str(len_ret/2 - jobs) \
                  + ' jobs and their settings will be ignored')

    return ret


def lower_dict_keys(dic):
    """ Turn all keys in a dictionary, or class derived from dictionary, to lowercase.
    Dictionaries are searched recursivly. """
    if isinstance(dic, dict):
        for key in dic:
            # Check if a key is lowercase; turn to lowercase if not
            try:
                key_lower = key.lower()
                if key != key_lower:
                    dic[key_lower] = dic[key]
                    del dic[key]
            except AttributeError:
                pass

            # Check if a value is a dictionary or a class derived from dictionary
            if isinstance(dic[key_lower], dict):
                dic[key_lower] = lower_dict_keys(dic[key_lower])

    return dic


def prep_qd(qd_list, path, arg):
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
        qd_list = list(QD_ligand_rotate.qd_opt(qd, qd_database, arg) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface
    if arg['qd_int']:
        print(QD_scripts.get_time() + 'calculating ligand distortion and inter-ligand interaction...')
        qd_list = list(QD_scripts.qd_int(qd) for qd in qd_list)

    # Calculate the interaction between ligands on the quantum dot surface upon removal of CdX2
    if arg['qd_dissociate']:
        # Extract the job type and input settings
        job1, s1, job2, s2 = get_job_settings(arg['qd_dissociate'], jobs=2)
        s1, s2 = lower_dict_keys(s1), lower_dict_keys(s2)

        # Start the BDE calculation
        print(QD_scripts.get_time() + 'calculating ligand dissociation energy...')
        for qd in qd_list:
            qd.properties.energy.BDE = QD_BDE.init_bde(qd, job1=job1, job2=job2, s1=s1, s2=s2)
            df = qd.properties.energy.BDE
            df.to_excel(join(path, qd.properties.name + '_BDE.xlsx'))

    # Write the new quantum dot results to the quantum dot database
    if arg['use_database']:
        if not arg['qd_opt']:
            for qd in qd_list:
                qd.properties.entry = True
        QD_database.write_database(qd_list, qd_database, path, mol_type='qd')

    return qd_list
