import os
import pandas as pd

import scm.plams.interfaces.molecule.rdkit as molkit
from rdkit import Chem
import QD_import_export as QD_inout


def read(mol_folder, database_name='ligand_database.xlsx'):
    """
    Open the database
    """
    database_path = os.path.join(mol_folder, database_name)
    if os.path.exists(database_path):
        database = pd.read_excel(database_path, sheet_name='Ligand')
    else:
        database = pd.DataFrame()

    return database


def compare(plams_mol, database):
    """
    Search the database for any ligand matches.
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
            QD_inout.set_prop(plams_mol_new, plams_mol.properties.name.replace('Ligand_', ''),
                              plams_mol.properties.source_folder)
            plams_mol = plams_mol_new
            pdb = True
        else:
            pdb = False
    else:
        match = False
        pdb = False

    return plams_mol, match, pdb


def write(ligand_list, database, database_name='ligand_database.xlsx'):
    """
    Write the new database entries to the database.
    """
    mol_folder = ligand_list[0].properties.source_folder
    database_entries = []
    for ligand in ligand_list:
        if ligand.properties.entry:
            prop = ligand.properties
            database_entries.append([prop.name,
                                     prop.formula,
                                     os.path.join(mol_folder, prop.name.split('@')[0]) + '.pdb',
                                     os.path.join(mol_folder, prop.name.split('@')[0]) + '.opt.pdb',
                                     prop.smiles])

    if database_entries:
        database_entries = list(zip(*database_entries))
        new = pd.DataFrame({'Ligand_name': database_entries[0],
                            'Ligand_formula': database_entries[1],
                            'Ligand_pdb': database_entries[2],
                            'Ligand_opt_pdb': database_entries[3],
                            'Ligand_SMILES': database_entries[4]})

        if not database.empty:
            new = database.append(new, ignore_index=True)

        database_path = os.path.join(mol_folder, database_name)
        new.to_excel(database_path, sheet_name='Ligand')
