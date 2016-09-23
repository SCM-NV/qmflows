# Default imports
from math import sqrt
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks import (Settings, templates)

# User Defined imports
from qmworks.packages.SCM import dftb
from qmworks.packages.orca import orca
from qmworks.packages import run


@attr('slow')
def test_opt_orca():
    """
    Test Orca input generation and run functions.
    """
    h2o = Molecule('test/test_files/h2o.xyz', 'xyz', charge=0, multiplicity=1)

    h2o_geometry = dftb(templates.geometry, h2o)

    s = Settings()
    # generic keyword "basis" must be present in the generic dictionary
    s.basis = "sto_dzp"
    # "specific" allows the user to apply specific keywords for a
    # package that are not in a generic dictionary
    # s.specific.adf.basis.core = "large"

    r = templates.singlepoint.overlay(s)
    h2o_singlepoint = orca(r, h2o_geometry.molecule)

    dipole = h2o_singlepoint.dipole

    final_result = run(dipole, n_processes=1)

    expected_dipole = [0.83203, 0.19521, -0.08393]

    diff = sqrt(sum((x - y) ** 2 for x, y in zip(final_result,
                                                 expected_dipole)))

    assert diff < 1e-8
