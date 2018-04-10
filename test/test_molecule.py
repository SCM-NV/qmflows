from nose.plugins.attrib import attr
from qmflows import molkit
from qmflows import packages


@attr('fast')
def test_SerMolecule():
    mol = molkit.from_smiles("c1ccccc1CC")
    registry = packages.registry()
    encoded_molecule = registry.deep_encode(mol)
    decoded_molecule = registry.deep_decode(encoded_molecule)
    assert len(mol) == len(decoded_molecule)
