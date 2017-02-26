
__all__ = ['add_prot_Hs', 'apply_reaction_smarts', 'apply_template', 'gen_coords_rdmol', 'modify_atom',
           'to_rdmol', 'from_rdmol', 'from_sequence', 'from_smiles', 'from_smarts', 'partition_protein',
           'write_molblock']

"""
@author: Lars Ridder
@description: A set of functions to manipulate molecules based on RDKit

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from plams import (Molecule, Bond, Atom)
import sys
import random
from warnings import warn


def from_rdmol(rdkit_mol, confid=-1):
    """
    Translate an RDKit molecule into a PLAMS molecule type.

    :parameter rdkit_mol: RDKit molecule
    :parameter int confid: conformer identifier from which to take coordinates
    :type rdkit_mol: rdkit.Chem.Mol
    :return: a PLAMS molecule
    :rtype: plams.Molecule

    """
    if isinstance(rdkit_mol, Molecule):
        return rdkit_mol
    # Create plams molecule
    plams_mol = Molecule()
    total_charge = 0
    Chem.Kekulize(rdkit_mol)
    conf = rdkit_mol.GetConformer(id=confid)
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
    Translate a PLAMS molecule into an RDKit molecule type.

    :parameter plams_mol: PLAMS molecule
    :type plams_mol: plams.Molecule
    :return: an RDKit molecule
    :rtype: rdkit.Chem.Mol

    """
    if isinstance(plams_mol, Chem.Mol):
        return plams_mol
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


def from_smiles(smiles, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates plams molecule(s) from a smiles strings.

    :parameter str smiles: A smiles string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization and ranking of comformations
                   The default value None results in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing similar/equivalent conformations
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    smiles = str(smiles.split()[0])
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdkit_mol.SetProp('smiles', smiles)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)

def from_smarts(smarts, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates plams molecule(s) from a smarts strings.
    This allows for example to define hydrogens explicitly.
    However it is less suitable for aromatic molecules (use from_smiles in that case).

    :parameter str smarts: A smarts string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization and ranking of comformations
                   The default value None results in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing similar/equivalent conformations
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    smiles = str(smarts.split()[0])
    mol = Chem.MolFromSmarts(smiles)
    Chem.SanitizeMol(mol)
    molecule = Chem.AddHs(mol)
    molecule.SetProp('smiles', smiles)
    return get_conformations(molecule, nconfs, name, forcefield, rms)

def get_conformations(rdkit_mol, nconfs=1, name=None, forcefield=None, rms=-1):
    """
    Generates 3D conformation(s) for and rdkit_mol

    :parameter rdkit_mol: RDKit molecule
    :type rdkit_mol: rdkit.Chem.Mol
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization and ranking of comformations
                   The default value None results in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing similar/equivalent conformations
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    def MMFFenergy(cid):
        ff = AllChem.MMFFGetMoleculeForceField(rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol), confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            warn("MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                 "\nNo geometry optimization was performed.")
            energy = 1e9
        return energy
    def UFFenergy(cid):
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            warn("MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                 "\nNo geometry optimization was performed.")
            energy = 1e9
        return energy

    if name:
        rdkit_mol.SetProp('name', name)
    cids = list(AllChem.EmbedMultipleConfs(rdkit_mol, nconfs, pruneRmsThresh=rms, randomSeed=1))
    if forcefield:
        optimize_molecule, energy = {
            'uff': [AllChem.UFFOptimizeMolecule, UFFenergy],
            'mmff': [AllChem.MMFFOptimizeMolecule, MMFFenergy],
            }[forcefield]
        for cid in cids:
            optimize_molecule(rdkit_mol, confId=cid)
        cids.sort(key=energy)
        if rms > 0:
            keep=[cids[0]]
            for cid in cids[1:]:
                for id in keep:
                    r = AllChem.AlignMol(rdkit_mol, rdkit_mol, cid, id)
                    if r < rms:
                        break
                else:
                    keep.append(cid)
            cids = keep
    if nconfs == 1:
        return from_rdmol(rdkit_mol)
    else:
        return [from_rdmol(rdkit_mol, cid) for cid in cids]

def from_sequence(sequence, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates plams molecule from a peptide sequence.
    Includes explicit hydrogens and 3D coordinates.

    :parameter str sequence: A peptide sequence, e.g. 'HAG'
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization and ranking of comformations
                   The default value None results in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing similar/equivalent conformations
    :return: A peptide molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    rdkit_mol = Chem.AddHs(Chem.MolFromSequence(sequence))
    rdkit_mol.SetProp('sequence', sequence)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)

def calc_rmsd(mol1, mol2):
    """
    Superimpose two molecules and calculate the root-mean-squared deviations of the atomic positions.
    The molecules should be identical, but the ordering of the atoms may differ.

    :param mol1: Molecule 1
    :param mol2: Molecule 2
    :return: The rmsd after superposition
    :rtype: float
    """
    rdkit_mol1 = to_rdmol(mol1)
    rdkit_mol2 = to_rdmol(mol2)
    try:
        return AllChem.GetBestRMS(rdkit_mol1, rdkit_mol2)
    except:
        return -999

def modify_atom(mol, idx, element):
    """
    Change atom "idx" in molecule "mol" to "element" and add or remove hydrogens accordingly

    :parameter mol: molecule to be modified
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :parameter int idx: index of the atom to be modified
    :parameter str element:
    :return: Molecule with new element and possibly added or removed hydrogens
    :rtype: plams.Molecule
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


def apply_template(mol, template):
    """
    Modifies bond orders in plams molecule according template smiles structure.

    :parameter mol: molecule to be modified
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :parameter str template: smiles string defining the correct chemical structure
    :return: Molecule with correct chemical structure and provided 3D coordinates
    :rtype: plams.Molecule
    """
    rdmol = to_rdmol(mol, sanitize=False)
    template_mol = Chem.AddHs(Chem.MolFromSmiles(template))
    newmol = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, rdmol)
    return from_rdmol(newmol)


def apply_reaction_smarts(mol, reaction_smarts, complete=False):
    """
    Applies reaction smirks and returns product.

    :parameter mol: molecule to be modified
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :parameter str reactions_smarts: Reactions smarts to be applied to molecule
    :parameter complete: Apply reaction until no further changes occur or given fraction of reaction centers have been modified
    :type complete: bool or float (value between 0 and 1)
    :return: (product molecule, list of unchanged atoms)
    :rtype: (plams.Molecule, list of int)
    """
    def react(reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        # if reaction doesn't apply, return the reactant
        if len(ps) == 0:
            return [(reactant, range(reactant.GetNumAtoms()))]
        full = len(ps)
        while complete: # when complete is True
            # apply reaction until no further changes
            r = random.randint(0,len(ps)-1)
            print(r)
            reactant = ps[r][0]
            ps = reaction.RunReactants([reactant])
            print('len:',len(ps))
            if len(ps) == 0 or len(ps)/full < (1-complete):
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
    """ Calculate 3D positions only for atoms without coordinates """
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
    Add hydrogens to protein molecules read from PDB.
    Makes sure that the hydrogens get the correct PDBResidue info.

    :param rdmol: An RDKit molecule containing a protein
    :type rdmol: rdkit.Chem.Mol
    :return: An RDKit molecule with explicit hydrogens added
    :rtype: rdkit.Chem.Mol
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


def partition_protein(rdmol):
    """
    Splits a protein molecule into capped amino acid fragments and caps.

    :param rdmol: A protein molecule
    :type rdmol: rdkit.Chem.Mol
    :return: list of fragments, list of caps
    """
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
