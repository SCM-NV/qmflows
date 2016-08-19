# Default imports
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks import (Settings, run)

# User Defined imports
from qmworks.components import (PES_scan, select_max)
from qmworks.packages.SCM import (dftb, adf)


@attr('slow')
def test_linear_ts():
    """
    compute a first approximation to the TS.
    """
    # Read the Molecule from file
    cnc = Molecule('C-N-C.mol', 'mol')

    # User define Settings
    settings = Settings()
    settings.functional = "pbe"
    settings.basis = "SZ"
    settings.specific.dftb.dftb.scc

    constraint1 = "dist 1 5"
    constraint2 = "dist 3 4"

    # scan input
    scan = {'constraint': [constraint1, constraint2],
            'surface': {'nsteps': 6, 'start': [2.3, 2.3],
                        'stepsize': [0.1, 0.1]}}

    # returns a set of results object containing the output of
    # each point in the scan
    lt = PES_scan([dftb, adf], settings, cnc, scan)

    # Gets the object presenting the molecule
    # with the maximum energy calculated from the scan
    apprTS = select_max(lt, "energy")

    # Run the TS optimization, using the default TS template
    ts = run(apprTS)

    expected_energy = -3.219708290363864

    assert abs(ts.energy - expected_energy) < 1.0e-6
    
