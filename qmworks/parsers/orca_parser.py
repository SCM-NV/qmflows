
from plams import (Atom, Molecule)
from pyparsing import (alphanums, Group, OneOrMore, Word)
from .parser import (floatNumber, skipLine, skipSupress,
                     string_array_to_molecule)
from .xyzParser import manyXYZ


def parse_molecule(file_name):
    """
    Parse The Cartesian coordinates from the output file.
    """
    header = "CARTESIAN COORDINATES (ANGSTROEM)"
    p1 = skipSupress(header) + skipLine * 2
    parse_atoms = Group(Word(alphanums) + floatNumber * 3)
    parse_mol = p1 + Group(OneOrMore(parse_atoms))

    parse_many_mol = OneOrMore(parse_mol)

    return string_array_to_molecule(parse_many_mol, file_name)


def parse_molecule_traj(file_traj):
    """
    Read Molecules from the *.traj file.
    """
    mols = manyXYZ(file_traj)
    # Last geometry corresponds to the optimized structure
    opt_mol = mols[-1]

    plams_mol = Molecule()
    for at in opt_mol:
        symb = at.symbol
        cs = at.xyz
        plams_mol.add_atom(Atom(symbol=symb, coords=tuple(cs)))

    return plams_mol


def parse_hessian(file_name):
    """ """
    pass
    
