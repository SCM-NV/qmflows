
from nose.plugins.attrib import attr
from noodles import gather
from qmworks import adf, run, Settings, templates
from qmworks.utils import initialize
import plams


@initialize
@attr('slow')
def test_ADFGeometry_Constraint():
    """
    Test "freeze" and "selected_atoms" keywords for constrained geometry
    optimizations.
    """
    an = plams.Molecule('test/test_files/an.xyz', 'xyz')
    # optimize only H atoms
    s = Settings()
    s.freeze = [1, 2, 3]
    result1 = adf(templates.geometry.overlay(s), an)
    geom1 = result1.molecule

    # optimize only H atoms
    s = Settings()
    s.selected_atoms = ['H']
    result2 = adf(templates.geometry.overlay(s), an)
    geom2 = result2.molecule

    r = run(gather(geom1, geom2), n_processes=1)

    assert str(r[0]) == str(r[1])


if __name__ == "__main__":
    plams.init()
    test_ADFGeometry_Constraint()
    plams.finish()
