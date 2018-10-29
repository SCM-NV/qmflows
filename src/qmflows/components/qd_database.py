__all__ = ['read_database', 'compare_database', 'write_database']

import os
import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem

from .qd_import_export import set_prop


def read_database(path, database_name='database.xlsx'):
    """
    Open the database.

    mol_folder <str>: The folder (including path) containing the database.
    database_name <str>: The name (including extension) of the database.

    return <pd.DataFrame>: A database of previous calculations.
    """
    path = os.path.join(path, database_name)
    if os.path.exists(path):
        database = pd.read_excel(path, sheet_name='Ligand')
    else:
        database = pd.DataFrame()

    return database


def compare_database(plams_mol, database):
    """
    Search the database for any ligand matches.

    plams_mol <plams.Molecule>: A plams molecule.
    database <pd.DataFrame>: A database of previous calculations.

    return <plams.Molecule>, <bool>, <bool>: The (imported) ligand, if a match was found between
        input ligand and the database, and if the .pdb file of this match actually exists.
    """
    mol_folder = plams_mol.properties.source_folder

    # If database usage is enabled: compare the ligand with previous entries.
    if not database.empty:
        database_mol = [Chem.MolFromSmiles(mol) for mol in database['Ligand_SMILES']]
        database_mol = [Chem.AddHs(mol) for mol in database_mol]
        matches = [molkit.to_rdmol(plams_mol).HasSubstructMatch(mol) for mol in database_mol]
    else:
        matches = [False]

    # Searches for matches between ligand and the database and Check if the ligand .pdb exists
    # Import the .pdb file
    if any(matches):
        index = matches.index(True)
        mol_path = os.path.join(mol_folder, str(database['Ligand_opt_pdb'][index]))
        match = True
        if os.path.exists(mol_path):
            plams_mol_new = molkit.readpdb(mol_path)
            prop_dict = {'mol_name': plams_mol.properties.name.replace('Ligand_', ''),
                         'folder_path': plams_mol.properties.source_folder}
            set_prop(plams_mol_new, prop_dict)
            plams_mol = plams_mol_new
            pdb = True
        else:
            pdb = False
    else:
        match = False
        pdb = False

    return plams_mol, match, pdb


def write_database(ligand_list, database, path, database_name='ligand_database.xlsx'):
    """
    Write the new database entries to the database.

    ligand_list <list>[<plams.Molecule>]: A list of ligands.
    database <pd.DataFrame>: A database of previous calculations.
    database_name <str>: The name (including extension) of the database.
    """
    mol_folder = ligand_list[0].properties.source_folder
    database_entries = []
    for ligand in ligand_list:
        if ligand.properties.entry:
            prop = ligand.properties
            database_entries.append(
                [prop.name,
                 prop.formula,
                 os.path.join(mol_folder, prop.name.split('@')[0]) + '.pdb',
                 os.path.join(mol_folder, prop.name.split('@')[0]) + '.opt.pdb',
                 prop.smiles,
                 prop.surface,
                 prop.volume,
                 prop.logp])

    if database_entries:
        database_entries = list(zip(*database_entries))
        new = pd.DataFrame({'Ligand_name': database_entries[0],
                            'Ligand_formula': database_entries[1],
                            'Ligand_pdb': database_entries[2],
                            'Ligand_opt_pdb': database_entries[3],
                            'Ligand_SMILES': database_entries[4],
                            'Ligand_surface': database_entries[5],
                            'Ligand_volume': database_entries[6],
                            'Ligand_logP': database_entries[7]})

        if not database.empty:
            new = database.append(new, ignore_index=True)

        path = os.path.join(path, database_name)
        new.to_excel(path, sheet_name='Ligand')


def write_database_qd(qd_list, path, database_name='qd_database.xlsx'):
    """
    Write the new database entries to the database.

    qd_list <list>[<plams.Molecule>]: A list of quantum_dots.
    database_name <str>: The name (including extension) of the database.
    """
    mol_folder = qd_list[0].properties.source_folder
    database_entries = []
    for qd in qd_list:
        if qd.properties:
            prop = qd.properties
            database_entries.append(
                [prop.name,
                 qd.get_formula(),
                 os.path.join(mol_folder, prop.name.split('@')[0]) + '.pdb',
                 os.path.join(mol_folder, prop.name.split('@')[0]) + '.opt.pdb',
                 prop.energy,
                 prop.int,
                 prop.strain])
    if database_entries:
        database_entries = list(zip(*database_entries))
        new = pd.DataFrame({'Quantum_dot_name': database_entries[0],
                            'Quantum_dot_formula': database_entries[1],
                            'Quantum_dot_pdb': database_entries[2],
                            'Quantum_dot_opt_pdb': database_entries[3],
                            'Quantum_dot_E': database_entries[4],
                            'Quantum_dot_Eint': database_entries[5],
                            'Quantum_dot_Estrain': database_entries[6]})

        path = os.path.join(path, database_name)
        if os.path.exists(path):
            database = pd.read_excel(path, sheet_name='Quantum_dot')
            database.append(new)
        new.to_excel(path, sheet_name='Quantum_dot')
