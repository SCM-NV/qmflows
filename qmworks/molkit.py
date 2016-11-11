
__all__ = ['apply_reaction_smarts', 'apply_template', 'gen_coords_rdmol', 'modify_atom',
           'to_rdmol', 'from_rdmol', 'from_sequence', 'from_smiles',
           'write_molblock']

"""
@author: Lars Ridder
@description: RDKit tools

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from plams import (Molecule, Bond, Atom)
import sys


def from_rdmol(rdkit_mol):
    """
    Translate an RDKit molecule into a PLAMS molecule type
    """
    plams_mol = Molecule()
    total_charge = 0
    Chem.Kekulize(rdkit_mol)
    conf = rdkit_mol.GetConformer()
    for atom in rdkit_mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        ch = atom.GetFormalCharge()
        plams_mol.add_atom(Atom(atom.GetAtomicNum(),
                                coords=(pos.x, pos.y, pos.z), charge=ch))
        total_charge += ch
    for bond in rdkit_mol.GetBonds():
        at1 = plams_mol.atoms[bond.GetBeginAtomIdx()]
        at2 = plams_mol.atoms[bond.GetEndAtomIdx()]
        plams_mol.add_bond(Bond(at1, at2, bond.GetBondTypeAsDouble()))
    plams_mol.charge = total_charge
    for propname in rdkit_mol.GetPropNames():
        plams_mol.properties[propname] = rdkit_mol.GetProp(propname)
    return plams_mol


def to_rdmol(plams_mol, sanitize=True):
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
        e.AddBond(a1, a2, Chem.BondType(bond.order))
    rdmol = e.GetMol()
    if sanitize:
        Chem.SanitizeMol(rdmol)
    conf = Chem.Conformer()
    for a in range(len(plams_mol.atoms)):
        atom = plams_mol.atoms[a]
        p = Geometry.Point3D(atom._getx(), atom._gety(), atom._getz())
        conf.SetAtomPosition(a, p)
    rdmol.AddConformer(conf)
    return rdmol


def from_smiles(smiles):
    """
    Generates plams molecule from a smiles strings.
    Includes explicit hydrogens and 3D coordinates
    """
    smiles = str(smiles.split()[0])
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    molecule.SetProp('smiles', smiles)
    AllChem.EmbedMolecule(molecule, randomSeed=1)
    AllChem.UFFOptimizeMolecule(molecule)
    return from_rdmol(molecule)


def from_sequence(sequence):
    """
    Generates plams molecule from a peptide sequence.
    Includes explicit hydrogens and 3D coordinates
    """
    molecule = Chem.MolFromSequence(sequence)
    AllChem.EmbedMolecule(molecule)
    AllChem.UFFOptimizeMolecule(molecule)
    return from_rdmol(Chem.AddHs(molecule, addCoords=True))


def modify_atom(mol, idx, element):
    """
    Change atom "idx" in molecule "mol" to "element"
    """
    rdmol = to_rdmol(mol)
    if rdmol.GetAtomWithIdx(idx).GetSymbol() == element:
        return mol
    else:
        e = Chem.EditableMol(rdmol)
        for neighbor in reversed(rdmol.GetAtomWithIdx(idx).GetNeighbors()):
            if neighbor.GetSymbol() == 'H':
                e.RemoveAtom(neighbor.GetIdx())
        e.ReplaceAtom(idx, Chem.Atom(element))
        newmol = e.GetMol()
        Chem.SanitizeMol(newmol)
        newmol = Chem.AddHs(newmol, addCoords=True)
        return from_rdmol(newmol)


def apply_template(plams_mol, template):
    """
    Modifies bond orders in plams molecule according template smiles structure
    """
    rdmol = to_rdmol(plams_mol, sanitize=False)
    template_mol = Chem.AddHs(Chem.MolFromSmiles(template))
    newmol = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, rdmol)
    return from_rdmol(newmol)


def apply_reaction_smarts(mol, reaction_smarts, complete=False):
    """
    Applies reaction smirks and returns product
    """
    def react(reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        # if reaction doesn't apply, return the reactant
        if len(ps) == 0:
            return [(reactant, range(reactant.GetNumAtoms()))]
        while complete: # when complete is True
            # apply reaction until no further changes
            reactant = ps[0][0]
            ps = reaction.RunReactants([reactant])
            print('len:',len(ps))
            if len(ps) == 0:
                ps = [[reactant]]
                break
        # add hydrogens and generate coordinates for new atoms
        products = []
        for p in ps[0]:
            Chem.SanitizeMol(p)
            q = Chem.AddHs(p)
            Chem.SanitizeMol(q)
            u = gen_coords_rdmol(q) # These are the atoms that have not changed
            products.append((q, u))
        return products

    if isinstance(mol, Molecule):
        mol = to_rdmol(mol)
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    # RDKit removes fragments that are disconnected from the reaction center
    # In order to keep these, the molecule is first split in separate fragments
    # and the results, including non-reacting parts, are re-combined afterwards
    frags = (Chem.GetMolFrags(mol, asMols=True))
    product = Chem.Mol()
    unchanged = [] # List of atoms that have not changed
    for frag in frags:
        for p, u in react(frag, reaction):
            unchanged += [product.GetNumAtoms() + i for i in u]
            product = Chem.CombineMols(product, p)
    # The molecule is returned together with a list of atom indices of the atoms that are identical to those
    # in the reactants. This list can be used in subsequent partial optimization of the molecule
    return (from_rdmol(product), unchanged)


def gen_coords(plamsmol):
    """ Calculate 3D positions for atoms without coordinates """
    rdmol = to_rdmol(plamsmol)
    unchanged = gen_coords_rdmol(rdmol)
    conf = rdmol.GetConformer()
    for a in range(len(plamsmol.atoms)):
        pos = conf.GetAtomPosition(a)
        atom = plamsmol.atoms[a]
        atom._setx(pos.x)
        atom._sety(pos.y)
        atom._setz(pos.z)
    return unchanged


def gen_coords_rdmol(rdmol):
    ref = rdmol.__copy__()
    conf = rdmol.GetConformer()
    coordDict = {}
    unchanged = []
    maps = []
    # Put known coordinates in coordDict
    for i in range(rdmol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if (-0.0001 < pos.x < 0.0001) and (-0.0001 < pos.y < 0.0001) and \
           (-0.0001 < pos.z < 0.0001):
            continue  # atom without coordinates
        coordDict[i] = pos
        unchanged.append(i)
        maps.append((i, i))
    # compute coordinates for new atoms, keeping known coordinates
    rms = 1
    rs = 1
    # repeat embedding and alignment until the rms of mapped atoms is sufficiently small
    while rms > 0.1:
        AllChem.EmbedMolecule(rdmol, coordMap=coordDict, randomSeed=rs,
                                    useBasicKnowledge=True)
        # align new molecule to original coordinates
        rms = AllChem.AlignMol(rdmol, ref, atomMap=maps)
        rs += 1
    return unchanged


def write_molblock(plams_mol, file=sys.stdout):
    file.write(Chem.MolToMolBlock(to_rdmol(plams_mol)))


def add_prot_Hs(rdmol):
    """
    Add hydrogens to molecules read from PDB
    Makes sure that the hydrogens get the correct PDBResidue info
    """
    retmol = Chem.AddHs(rdmol, addCoords=True)
    for atom in retmol.GetAtoms():
        if atom.GetPDBResidueInfo() is None and atom.GetSymbol() == 'H':
            bond = atom.GetBonds()[0]
            if bond.GetBeginAtom().GetIdx() == atom.GetIdx:
                connected_atom = bond.GetEndAtom()
            else:
                connected_atom = bond.GetBeginAtom()
            try:
                ResInfo = connected_atom.GetPDBResidueInfo()
                atom.SetMonomerInfo(ResInfo)
            except:
                print('Hydrogen annotation failed:', connected_atom.GetIdx(),
                      atom.GetIdx())
    return retmol


def add_fragment(rwmol, frag, rwmol_atom_idx=None, frag_atom_idx=None,
                 bond_order=None):
    molconf = rwmol.GetConformer()
    fragconf = frag.GetConformer()
    new_indices = []
    for a in frag.GetAtoms():
        new_index = rwmol.AddAtom(a)
        new_indices.append(new_index)
        molconf.SetAtomPosition(new_index, fragconf.GetAtomPosition(a.GetIdx()))
    for b in frag.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        rwmol.AddBond(new_indices[ba], new_indices[ea], b.GetBondType())
    if bond_order:
        rwmol.AddBond(rwmol_atom_idx, new_indices[frag_atom_idx],
                      Chem.BondType.values[bond_order])
        rwmol.GetAtomWithIdx(new_indices[frag_atom_idx]).SetNumRadicalElectrons(0)


def get_fragment(mol, indices, incl_expl_Hs=True, neutralize=True):
    molconf = mol.GetConformer()
    fragment = Chem.RWMol(Chem.Mol())
    fragconf = Chem.Conformer()
    # Put atoms in fragment
    for i in indices:
        atom = mol.GetAtomWithIdx(i)
        new_index = fragment.AddAtom(atom)
        pos = molconf.GetAtomPosition(i)
        fragconf.SetAtomPosition(new_index, pos)
    # Put bonds in fragment
    for b in mol.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        if ba in indices and ea in indices:
            fragment.AddBond(indices.index(ba), indices.index(ea),
                             b.GetBondType())
            continue
        if not incl_expl_Hs:
            continue
        if ba in indices and mol.GetAtomWithIdx(ea).GetSymbol() == 'H':
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ea))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ea))
            fragment.AddBond(indices.index(ba), hi, Chem.BondType.SINGLE)
            continue
        if ea in indices and mol.GetAtomWithIdx(ba).GetSymbol() == 'H':
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ba))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ba))
            fragment.AddBond(indices.index(ea), hi, Chem.BondType.SINGLE)
    ret_frag = fragment.GetMol()
    Chem.SanitizeMol(ret_frag)
    if neutralize:
        for atom in ret_frag.GetAtoms():
            nrad = atom.GetNumRadicalElectrons()
            if nrad > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + nrad)
                atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(ret_frag)
    ret_frag.AddConformer(fragconf)
    return ret_frag


def partition_protein(rdmol, cap=None):
    caps = []
    em = Chem.RWMol(rdmol)
    # Split peptide bonds
    pept_bond = Chem.MolFromSmarts('[C;X4;H1,H2][CX3](=O)[NX3][C;X4;H1,H2][CX3](=O)')
    for match in rdmol.GetSubstructMatches(pept_bond):
        cap = get_fragment(rdmol, match[0:5])
        cap = add_prot_Hs(cap)
        caps.append(cap)
        cap_o_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4][CX3]=O'))
        cap_o = get_fragment(cap, cap_o_ind, neutralize=False)
        cap_n_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('O=[CX3][NX3][C;X4]'))[2:]
        cap_n = get_fragment(cap, cap_n_ind, neutralize=False)
        em.RemoveBond(match[1], match[3])
        add_fragment(em, cap_o, match[3], 1, 1)
        add_fragment(em, cap_n, match[1], 0, 1)
    # Split disulfide bonds
    ss_bond = Chem.MolFromSmarts('[C;X4;H1,H2]SS[C;X4;H1,H2]')
    for match in rdmol.GetSubstructMatches(ss_bond):
        cap = get_fragment(rdmol, match[0:5])
        cap = add_prot_Hs(cap)
        caps.append(cap)
        cap_s_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4]SS[C;X4]'))
        cap_s1 = get_fragment(cap, cap_s_ind[0:2], neutralize=False)
        cap_s2 = get_fragment(cap, cap_s_ind[2:4], neutralize=False)
        em.RemoveBond(match[1], match[2])
        add_fragment(em, cap_s1, match[2], 1, 1)
        add_fragment(em, cap_s2, match[1], 0, 1)
    Chem.SanitizeMol(em)
    frags = Chem.GetMolFrags(em.GetMol(), asMols=True)
    return frags, caps
