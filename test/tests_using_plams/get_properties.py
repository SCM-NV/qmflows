from qmworks import (templates, run, rdkitTools, dftb)

@attr('slow')
def test_dftb_props():
    """
    Get properties from DFTB freq calc
    """

    mol = rdkitTools.smiles2rdkit('F[H]')
    result = run(dftb(templates.freq, mol, job_name='dftb_FH'))
    expected_energy = -4.76
    assert abs(result.energy - expected_energy) < 0.01
    assert len(result.dipole) == 3
    expected_frequency = 3293.85
    assert abs(result.frequencies - expected_frequency) < 0.01
    assert len(result.charges) == 2
    assert abs(result.charges[0] + result.charges[1]) < 1e-6