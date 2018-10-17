from qmflows import packages
import pytest
import scm.plams.interfaces.molecule.rdkit as molkit


@pytest.mark.slow
def test_SerMolecule():
    mol = molkit.from_smiles("c1ccccc1CC")
    registry = packages.registry()
    encoded_molecule = registry.deep_encode(mol)
    decoded_molecule = registry.deep_decode(encoded_molecule)
    assert len(mol) == len(decoded_molecule)
