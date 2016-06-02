"""
@author: Lars Ridder
@description: RDKit tools

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, rdMolTransforms
from plams import Molecule, Bond, Atom
from io import StringIO
import sys

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
        ch = atom.GetFormalCharge()
        plams_mol.add_atom(Atom(atom.GetAtomicNum(),coords = (pos.x,pos.y,pos.z), charge=ch))
        total_charge += ch
    for bond in rdkit_mol.GetBonds():
        at1 = plams_mol.atoms[bond.GetBeginAtomIdx()]
        at2 = plams_mol.atoms[bond.GetEndAtomIdx()]
        plams_mol.add_bond(Bond(at1,at2,bond.GetBondTypeAsDouble()))
    plams_mol.charge = total_charge
    return plams_mol

def plams2rdkit(plams_mol, sanitize=True):
    """
    Translate a PLAMS molecule into an RDKit molecule type
    """
    # Create rdkit molecule
    e = Chem.EditableMol(Chem.Mol())
    for atom in plams_mol.atoms:
        a = Chem.Atom(atom.atnum)
        ch = atom.properties.charge
        if isinstance(ch, int):
            a.SetFormalCharge(ch)
        e.AddAtom(a)
    for bond in plams_mol.bonds:
        a1 = plams_mol.atoms.index(bond.atom1)
        a2 = plams_mol.atoms.index(bond.atom2)
        e.AddBond(a1,a2,Chem.BondType(bond.order))
    rdmol = e.GetMol()
    if sanitize:
        Chem.SanitizeMol(rdmol)
    conf = Chem.Conformer()
    for a in range(len(plams_mol.atoms)):
        atom = plams_mol.atoms[a]
        p = Geometry.Point3D(atom._getx(), atom._gety(), atom._getz())
        conf.SetAtomPosition(a,p)
    rdmol.AddConformer(conf)    
    return rdmol
       
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

def apply_template(plams_mol, template):
    """
    Modifies bond orders in plams molecule according template smiles structure
    """
    rdmol = plams2rdkit(plams_mol, sanitize=False)
    template_mol = Chem.AddHs(Chem.MolFromSmiles(template))
    newmol = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol,rdmol)
    return rdkit2plams(newmol)

def apply_smirks(plams_mol, reaction_smirks):
    """
    Applies reaction smirks and returns list of products
    """
    def react(reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        products = []
        for product in ps:
            frags = (Chem.GetMolFrags(product[0], asMols=True))
            for p in frags:
                q = Chem.AddHs(p)
                Chem.SanitizeMol(q)
                #gen_coords(q)
                products.append(q)
        return products

    rdmol = plams2rdkit(plams_mol)
    query = Chem.MolFromSmarts(reaction_smirks.split('>>')[0])
    reaction = AllChem.ReactionFromSmarts(reaction_smirks)
    products = react(rdmol, reaction)
    product_list = [rdkit2plams(p) for p in products]
    return product_list

def gen_coords(plams_mol):
    """ Calculate 3D positions for atoms without coordinates """
    rdmol = plams2rdkit(plams_mol)
    conf = rdmol.GetConformer(0)
    coordDict = {}
    # Put known coordinates in coordDict
    for i in range(rdmol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if (-0.0001 < pos.x < 0.0001) and (-0.0001 < pos.y < 0.0001) and (-0.0001 < pos.z < 0.0001):
            continue # atom without coordinates
        coordDict[i] = Geometry.Point3D(pos.x, pos.y, pos.z)
    # compute coordinates for new atoms, keeping known coordinates
    seed = 1
    i = -1
    while i == -1:
        i = AllChem.EmbedMolecule(rdmol,
                                  ignoreSmoothingFailures=True,
                                  randomSeed=seed,
                                  useBasicKnowledge=True,
                                  coordMap=coordDict,
                                  forceTol=0.01)
        seed += 1
    return rdkit2plams(rdmol)

def write_molblock(plams_mol, file=sys.stdout):
    file.write(Chem.MolToMolBlock(plams2rdkit(plams_mol)))
