
__all__ = ['add_Hs', 'apply_reaction_smarts', 'apply_template',
           'gen_coords_rdmol', 'get_backbone_atoms', 'modify_atom',
           'to_rdmol', 'from_rdmol', 'from_sequence', 'from_smiles',
           'from_smarts', 'partition_protein', 'readpdb', 'writepdb']

"""
@author: Lars Ridder
@description: A set of functions to manipulate molecules based on RDKit

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from scm.plams import (Molecule, Bond, Atom)
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
    try:
        Chem.Kekulize(rdkit_mol)
    except:
        pass
    conf = rdkit_mol.GetConformer(id=confid)
    for rd_atom in rdkit_mol.GetAtoms():
        pos = conf.GetAtomPosition(rd_atom.GetIdx())
        ch = rd_atom.GetFormalCharge()
        pl_atom = Atom(
            rd_atom.GetAtomicNum(), coords=(pos.x, pos.y, pos.z), charge=ch)
        if rd_atom.GetPDBResidueInfo():
            pl_atom.properties.pdb_info = get_PDBResidueInfo(rd_atom)
        plams_mol.add_atom(pl_atom)
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
        if 'charge' in atom.properties:
            a.SetFormalCharge(atom.properties.charge)
        if 'pdb_info' in atom.properties:
            set_PDBresidueInfo(a, atom.properties.pdb_info)
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


pdb_residue_info_items = [
    'AltLoc', 'ChainId', 'InsertionCode', 'IsHeteroAtom', 'Name', 'Occupancy',
    'ResidueName', 'ResidueNumber', 'SecondaryStructure', 'SegmentNumber',
    'SerialNumber', 'TempFactor']
# 'MonomerType' was excluded because it is an rdkit type that cannot easilty be serialized


def get_PDBResidueInfo(rdkit_atom):
    pdb_info = {}
    for item in pdb_residue_info_items:
        get_function = 'Get' + item
        pdb_info[item] = rdkit_atom.GetPDBResidueInfo().__getattribute__(get_function)()
    return pdb_info


def set_PDBresidueInfo(rdkit_atom, pdb_info):
    atom_pdb_residue_info = Chem.AtomPDBResidueInfo()
    for item, value in pdb_info.items():
        set_function = 'Set' + item
        atom_pdb_residue_info.__getattribute__(set_function)(value)
    rdkit_atom.SetMonomerInfo(atom_pdb_residue_info)


def from_smiles(smiles, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates plams molecule(s) from a smiles strings.

    :parameter str smiles: A smiles string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization
        and ranking of comformations. The default value None results in skipping of the
        geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for
        removing similar/equivalent conformations
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
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
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
    Generates 3D conformation(s) for an rdkit_mol

    :parameter rdkit_mol: RDKit molecule
    :type rdkit_mol: rdkit.Chem.Mol
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
    optimization and ranking of comformations. The default value None results
    in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing
    similar/equivalent conformations.
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    def MMFFenergy(cid):
        ff = AllChem.MMFFGetMoleculeForceField(
            rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol), confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = "MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                  "\nNo geometry optimization was performed."
            warn(msg)
            energy = 1e9
        return energy

    def UFFenergy(cid):
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = "MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                  "\nNo geometry optimization was performed."
            warn(msg)
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
            keep = [cids[0]]
            for cid in cids[1:]:
                for idx in keep:
                    try:
                        r = AllChem.AlignMol(rdkit_mol, rdkit_mol, cid, idx)
                    except:
                        r = rms + 1
                        message = "Alignment failed in multiple conformation generation: "
                        message += Chem.MolToSmiles(rdkit_mol)
                        message += "\nAssuming different conformations."
                        warn(message)
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
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :return: A peptide molecule with hydrogens and 3D coordinates
        or a list of molecules if nconfs > 1
    :rtype: plams.Molecule or list of plams Molecules
    """
    rdkit_mol = Chem.AddHs(Chem.MolFromSequence(sequence))
    rdkit_mol.SetProp('sequence', sequence)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)


def calc_rmsd(mol1, mol2):
    """
    Superimpose two molecules and calculate the root-mean-squared deviations of
    the atomic positions.
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
        for neighbor in reversed(rdmol.GetAtomWithIdx(idx - 1).GetNeighbors()):
            if neighbor.GetSymbol() == 'H':
                e.RemoveAtom(neighbor.GetIdx())
        e.ReplaceAtom(idx - 1, Chem.Atom(element))
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


def apply_reaction_smarts(
        mol, reaction_smarts, complete=False, forcefield=None, return_rdmol=False):
    """
    Applies reaction smirks and returns product.
    If returned as a plams molecule, plams.Molecule.properties.orig_atoms
    is a list of indices of atoms that have not been changed
    (which can for example be used partially optimize new atoms only with the freeze keyword)

    :parameter mol: molecule to be modified
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :parameter str reactions_smarts: Reactions smarts to be applied to molecule
    :parameter complete: Apply reaction until no further changes occur or given
        fraction of reaction centers have been modified
    :type complete: bool or float (value between 0 and 1)
    :parameter forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization of product structures.
    :type forcefield: str
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
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
        while complete:  # when complete is True
            # apply reaction until no further changes
            r = random.randint(0, len(ps) - 1)
            reactant = ps[r][0]
            ps = reaction.RunReactants([reactant])
            if len(ps) == 0 or len(ps) / full < (1 - complete):
                ps = [[reactant]]
                break
        # add hydrogens and generate coordinates for new atoms
        products = []
        for p in ps[0]:
            Chem.SanitizeMol(p)
            q = Chem.AddHs(p)
            Chem.SanitizeMol(q)
            u = gen_coords_rdmol(q)  # These are the atoms that have not changed
            products.append((q, u))
        return products

    mol = to_rdmol(mol)
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    # RDKit removes fragments that are disconnected from the reaction center
    # In order to keep these, the molecule is first split in separate fragments
    # and the results, including non-reacting parts, are re-combined afterwards
    frags = (Chem.GetMolFrags(mol, asMols=True))
    product = Chem.Mol()
    unchanged = []  # List of atoms that have not changed
    for frag in frags:
        for p, u in react(frag, reaction):
            unchanged += [product.GetNumAtoms() + i for i in u]
            product = Chem.CombineMols(product, p)
    if forcefield:
        optimize_coordinates(product, forcefield, fixed=unchanged)
    # The molecule is returned together with a list of atom indices of the atoms
    # that are identical to those
    # in the reactants. This list can be used in subsequent partial optimization of the molecule
    if not return_rdmol:
        product = from_rdmol(product)
        product.properties.orig_atoms = [a + 1 for a in unchanged]
    return product


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
    return [a + 1 for a in unchanged]


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
    if rdmol.GetNumAtoms() > len(maps):
        while rms > 0.1:
            AllChem.EmbedMolecule(rdmol, coordMap=coordDict, randomSeed=rs,
                                  useBasicKnowledge=True)
            # align new molecule to original coordinates
            rms = AllChem.AlignMol(rdmol, ref, atomMap=maps)
            rs += 1
    return unchanged


def optimize_coordinates(rdkit_mol, forcefield, fixed=[]):
    def MMFFminimize():
        ff = AllChem.MMFFGetMoleculeForceField(
            rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol))
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("MMFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))

    def UFFminimize():
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, ignoreInterfragInteractions=True)
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("UFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))
    optimize_molecule = {
        'uff': UFFminimize,
        'mmff': MMFFminimize}[forcefield]
    Chem.SanitizeMol(rdkit_mol)
    optimize_molecule()
    return


def write_molblock(plams_mol, file=sys.stdout):
    file.write(Chem.MolToMolBlock(to_rdmol(plams_mol)))


def readpdb(pdb_file, removeHs=False, return_rdmol=False):
    """
    Generate a molecule from a PDB file

    :param pdb_file: The PDB file to read
    :type pdb_file: str or file
    :param bool removeHs: Hydrogens are removed if Trur
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: The molecule
    :rtype: plams.Molecule or rdkit.Chem.Mol
    """
    if isinstance(pdb_file, str):
        pdb_file = open(pdb_file, 'r')
    pdb_mol = Chem.MolFromPDBBlock(pdb_file.read(), removeHs=removeHs)
    return pdb_mol if return_rdmol else from_rdmol(pdb_mol)


def writepdb(mol, pdb_file=sys.stdout):
    """
    Write a PDB file from a molecule

    :parameter mol: molecule to be exported to PDB
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :param pdb_file: The PDB file to write to, or a filename
    :type pdb_file: str or file
    """
    mol = to_rdmol(mol)
    if isinstance(pdb_file, str):
        pdb_file = open(pdb_file, 'w')
    pdb_file.write(Chem.MolToPDBBlock(mol))


def add_Hs(mol, forcefield=None, return_rdmol=False):
    """
    Add hydrogens to protein molecules read from PDB.
    Makes sure that the hydrogens get the correct PDBResidue info.

    :param mol: Molecule to be protonated
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :param str forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization on new atoms.
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: A molecule with explicit hydrogens added
    :rtype: plams.Molecule or rdkit.Chem.Mol
    """
    mol = to_rdmol(mol)
    retmol = Chem.AddHs(mol)
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
                atomname = 'H' + atom.GetPDBResidueInfo().GetName()[1:]
                atom.GetPDBResidueInfo().SetName(atomname)
            except:
                pass
    unchanged = gen_coords_rdmol(retmol)
    if forcefield:
        optimize_coordinates(retmol, forcefield, fixed=unchanged)
    return retmol if return_rdmol else from_rdmol(retmol)


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


def partition_protein(mol, residue_bonds=None, split_heteroatoms=True, return_rdmol=False):
    """
    Splits a protein molecule into capped amino acid fragments and caps.

    :param mol: A protein molecule
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :param tuple residue_bonds: a tuple of pairs of residue number indicating which
        peptide bonds to split. If none, split all peptide bonds.
    :param bool split_heteroatoms: if True, all bonds between a heteroatom and
        a non-heteroatom across residues are removed
    :return: list of fragments, list of caps
    """
    mol = to_rdmol(mol)
    caps = []
    em = Chem.RWMol(mol)
    if split_heteroatoms:
        for bond in mol.GetBonds():
            resinfa = bond.GetBeginAtom().GetPDBResidueInfo()
            resinfb = bond.GetEndAtom().GetPDBResidueInfo()
            if resinfa.GetIsHeteroAtom() is not resinfb.GetIsHeteroAtom():
                if resinfa.GetResidueNumber() != resinfb.GetResidueNumber():
                    em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    # Split peptide bonds
    pept_bond = Chem.MolFromSmarts('[C;X4;H1,H2][CX3](=O)[NX3][C;X4;H1,H2][CX3](=O)')
    for match in mol.GetSubstructMatches(pept_bond):
        if residue_bonds:
            resa = mol.GetAtomWithIdx(match[1]).GetPDBResidueInfo().GetResidueNumber()
            resb = mol.GetAtomWithIdx(match[3]).GetPDBResidueInfo().GetResidueNumber()
            if (resa, resb) not in residue_bonds and (resb, resa) not in residue_bonds:
                continue
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_o_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4][CX3]=O'))
        cap_o = get_fragment(cap, cap_o_ind, neutralize=False)
        cap_n_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('O=[CX3][NX3][C;X4]'))[2:]
        cap_n = get_fragment(cap, cap_n_ind, neutralize=False)
        em.RemoveBond(match[1], match[3])
        add_fragment(em, cap_o, match[3], 1, 1)
        add_fragment(em, cap_n, match[1], 0, 1)
    # Split disulfide bonds
    ss_bond = Chem.MolFromSmarts('[C;X4;H1,H2]SS[C;X4;H1,H2]')
    for match in mol.GetSubstructMatches(ss_bond):
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_s_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4]SS[C;X4]'))
        cap_s1 = get_fragment(cap, cap_s_ind[0:2], neutralize=False)
        cap_s2 = get_fragment(cap, cap_s_ind[2:4], neutralize=False)
        em.RemoveBond(match[1], match[2])
        add_fragment(em, cap_s1, match[2], 1, 1)
        add_fragment(em, cap_s2, match[1], 0, 1)
    frags = Chem.GetMolFrags(em.GetMol(), asMols=True, sanitizeFrags=False)
    if not return_rdmol:
        frags = [from_rdmol(frag) for frag in frags]
    return frags, caps


def charge_AAs(mol, return_rdmol=False):
    ionizations = {
        'ARG_NH2': 1,
        'LYS_NZ': 1,
        'GLU_OE2': -1,
        'ASP_OD2': -1}
    mol = to_rdmol(mol)
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        res_atom = resinfo.GetResidueName() + '_' + resinfo.GetName().strip()
        try:
            atom.SetFormalCharge(ionizations[res_atom])
            Chem.SanitizeMol(mol)
        except KeyError:
            pass
        Chem.SanitizeMol(mol)
    return mol if return_rdmol else from_rdmol(mol)


def get_backbone_atoms(mol):
    """
    Return a list of atom indices corresponding to the backbone atoms in a peptide molecule.
    This function assumes PDB information in properties.pdb_info of each atom, which is the case
    if the molecule is generated with the "readpdb" or "from_sequence" functions.

    :parameter mol: a peptide molecule
    :type mol: plams.Molecule or rdkit.Chem.Mol
    :return: a list of atom indices
    :rtype: list
    """
    mol = from_rdmol(mol)
    backbone = ['N', 'CA', 'C', 'O']
    return [a for a in range(1, len(mol) + 1)
            if str(mol[a].properties.pdb_info.Name).strip() in backbone]
