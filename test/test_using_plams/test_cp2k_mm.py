import pytest
from assertionlib import assertion
from scm.plams import Molecule, add_to_instance, Units

from qmflows import Settings, run, cp2k_mm, singlepoint, geometry
from qmflows.test_utils import delete_output, get_mm_settings, PATH

MOL = Molecule(PATH / 'Cd68Se55_26COO.xyz')
MOL.guess_bonds([at for at in MOL if at.symbol not in {'Cd', 'Se'}])

SETTINGS: Settings = get_mm_settings()


def get_energy(self, index: int = -1, unit: str = 'Hartree') -> float:
    """Return the energy of the last occurence of ``'ENERGY| Total FORCE_EVAL'`` in the output."""
    energy_str = self.grep_output('ENERGY| Total FORCE_EVAL')[index]
    energy = float(energy_str.rsplit(maxsplit=1)[1])
    return Units.convert(energy, 'Hartree', unit)


@pytest.mark.slow
@delete_output
def test_singlepoint() -> None:
    """Test CP2K singlepoint calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += singlepoint.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Yes, this is a small hack as neither energy nor get_energy() seems to work
    plams_results = result.results
    add_to_instance(plams_results)(get_energy)
    setattr(result, 'get_energy', plams_results.get_energy)

    ref = -32.390062964050124
    energy = result.get_energy()
    assertion.isclose(energy, ref, rel_tol=10**-4)


@pytest.mark.slow
@delete_output(delete_workdir=False)
def test_geometry() -> None:
    """Test CP2K geometry optimization calculations with the :class:`CP2K_MM` class."""
    s = SETTINGS.copy()
    s.specific.cp2k += geometry.specific.cp2k_mm

    job = cp2k_mm(settings=s, mol=MOL)
    result = run(job, path=PATH)
    assertion.eq(result.status, 'successful')

    # Yes, this is a small hack as neither energy nor get_energy() seems to work
    plams_results = result.results
    add_to_instance(plams_results)(get_energy)
    setattr(result, 'get_energy', plams_results.get_energy)

    ref = -32.390062964050124
    energy = result.get_energy()
    assertion.isclose(energy, ref, rel_tol=10**-4)
