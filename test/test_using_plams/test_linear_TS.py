# Default imports
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks import (Settings, run)

# User Defined imports
from qmworks.components import (PES, select_max, Distance)
from qmworks.packages.SCM import (dftb, adf)


@attr('slow')
def test_linear_ts():
    """
    compute a first approximation to the TS.
    """
    # Read the Molecule from file
    cnc = Molecule('test/test_files/C-N-C.mol', 'mol')

    # User define Settings
    settings = Settings()
    settings.functional = "pbe"
    settings.basis = "SZ"
    settings.specific.dftb.dftb.scc

    constraint1 = Distance(0, 4)
    constraint2 = Distance(2, 3)

    # scan input
    pes = PES(cnc, constraints=[constraint1, constraint2],
              offset=[2.3, 2.3], get_current_values=False, nsteps=2, stepsize=[0.1, 0.1])

    # returns a set of results object containing the output of
    # each point in the scan
    lt = pes.scan([dftb, adf], settings)

    # Gets the object presenting the molecule
    # with the maximum energy calculated from the scan
    apprTS = select_max(lt, "energy")

    # Run the TS optimization, using the default TS template
    ts = run(apprTS)

    expected_energy = -3.219708290363864

    assert abs(ts.energy - expected_energy) < 0.02
