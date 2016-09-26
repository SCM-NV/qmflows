
from pyparsing import (alphanums, Group, OneOrMore, Word)
from .parser import (floatNumber, skipLine, skipSupress, string_array_to_molecule)


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
