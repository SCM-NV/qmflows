"""
@author: Lars Ridder
@description: RDKit tools

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from plams import Molecule,Bond,Atom
from io import StringIO

def rdkit2plams(rdkit_mol):
    """
    Translate an RDKit molecule into a PLAMS molecule type
    """
    plams_mol = Molecule()
    total_charge = 0
    Chem.Kekulize(rdkit_mol)
    conf=rdkit_mol.GetConformer()
    for atom in rdkit_mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        plams_mol.add_atom(Atom(atom.GetAtomicNum(),coords = (pos.x,pos.y,pos.z)))
        total_charge += atom.GetFormalCharge()
    for bond in rdkit_mol.GetBonds():
        at1 = plams_mol.atoms[bond.GetBeginAtomIdx()]
        at2 = plams_mol.atoms[bond.GetEndAtomIdx()]
        plams_mol.add_bond(Bond(at1,at2,bond.GetBondTypeAsDouble()))
    plams_mol.charge = total_charge
    return plams_mol

def plams2rdkit(plams_mol):
    """
    Translate a PLAMS molecule into an RDKit molecule type
    """
    molblock = StringIO()
    plams_mol.writemol(molblock)
    rdkit_mol=Chem.MolFromMolBlock(molblock.getvalue(),removeHs=False)
    return rdkit_mol
        
def smiles2plams(smiles):
    """
    Generates plams molecules from a smiles strings.
    Includes explicit hydrogens and 3D coordinates
    """
    smiles = str(smiles.split()[0])
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed = 1)
    AllChem.UFFOptimizeMolecule(molecule)
    return rdkit2plams(molecule)
#    return rdkit2plams(Chem.AddHs(molecule,addCoords=True))

def sequence2plams(sequence):
    """
    Generates plams molecules from a peptide sequence.
    Includes explicit hydrogens and 3D coordinates
    """
    molecule = Chem.MolFromSequence(sequence)
    AllChem.EmbedMolecule(molecule)
    AllChem.UFFOptimizeMolecule(molecule)
    print(Chem.MolToMolBlock(molecule))
    return rdkit2plams(Chem.AddHs(molecule,addCoords=True))

def modify_atom(mol, idx, element):
    """
    Change atom "idx" in molecule "mol" to "element"
    """
    rdmol = plams2rdkit(mol)
    if rdmol.GetAtomWithIdx(idx).GetSymbol() == element:
        return mol
    else:
        e = Chem.EditableMol(rdmol)
        for neighbor in reversed(rdmol.GetAtomWithIdx(idx).GetNeighbors()):
            if neighbor.GetSymbol() == 'H':
                e.RemoveAtom(neighbor.GetIdx())
        e.ReplaceAtom(idx,Chem.Atom(element))
        newmol=e.GetMol()
        Chem.SanitizeMol(newmol)
        newmol=Chem.AddHs(newmol,addCoords=True)
        return rdkit2plams(newmol)

