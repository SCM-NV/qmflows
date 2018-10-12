import pandas as pd
import os

from qmflows import molkit
from rdkit import Chem


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


def compare(plams_mol, mol_name, mol_folder, database):
    """
    Search the database for any ligand matches.
    """
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
            plams_mol = molkit.readpdb(mol_path)
            pdb = True
        else:
            pdb = False
    else:
        match = False
        pdb = False

    return plams_mol, match, pdb


def entry(plams_mol, mol_name, mol_smiles, opt):
    """
    Create a new entry for the database.
    """
    a = mol_name
    b = plams_mol.get_formula()
    c = a + '.pdb'
    if opt:
        d = a + '.opt.pdb'
    else:
        d = 'ligand_optimization_disabled'
    e = mol_smiles

    return [a, b, c, d, e]


def write(database_entries, database, mol_folder, database_name='ligand_database.xlsx'):
    """
    Write the new database entries to the database.
    """
    database_entries = [item for item in database_entries if item]
    if database_entries:
        database_entries = list(zip(*database_entries))
        database_path = os.path.join(mol_folder, database_name)
        new = pd.DataFrame({'Ligand_name' : database_entries[0],
                            'Ligand_formula' : database_entries[1],
                            'Ligand_pdb' : database_entries[2],
                            'Ligand_opt_pdb' : database_entries[3],
                            'Ligand_SMILES' : database_entries[4]})
    
        if not database.empty:
            new = database.append(new, ignore_index=True)

        new.to_excel(database_path, sheet_name='Ligand')
