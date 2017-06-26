__all__ = ['example_FDE_fragments']

# Default imports
from scm.plams import Molecule
from qmworks import (Settings, run, templates)
from qmworks.packages.SCM import adf
from qmworks.components import Fragment, adf_fragmentsjob
from noodles import gather
import io


def example_FDE_fragments():
    # For the purpose of the example, define xyz files here:
    xyz1 = io.StringIO('''3

    O         0.00000000000000     -2.29819386240000      1.63037963360000
    H        -0.76925379540000     -2.28223123190000      2.22684542850000
    H         0.76925379540000     -2.28223123190000      2.22684542850000''')

    xyz2 = io.StringIO('''3

    O         0.00000000000000      2.29819386240000      1.63037963360000
    H        -0.76925379540000      2.28223123190000      2.22684542850000
    H         0.76925379540000      2.28223123190000      2.22684542850000''')

    xyz3 = io.StringIO('''3

    O         0.00000000000000      0.00000000000000     -0.26192472620000
    H         0.00000000000000      0.77162768440000      0.34261631290000
    H         0.00000000000000     -0.77162768440000      0.34261631290000''')

    # Read the Molecule from file
    m_h2o_1 = Molecule()
    m_h2o_1.readxyz(xyz1, 1)
    m_h2o_2 = Molecule()
    m_h2o_2.readxyz(xyz2, 1)
    m_mol = Molecule()
    m_mol.readxyz(xyz3, 1)

    settings = Settings()
    settings.basis = 'SZ'
    settings.specific.adf.nosymfit = ''

    # Prepare first water fragment
    r_h2o_1 = adf(templates.singlepoint.overlay(settings), m_h2o_1, job_name="h2o_1")

    # Prepare second water fragment
    r_h2o_2 = adf(templates.singlepoint.overlay(settings), m_h2o_2, job_name="h2o_2")

    frags = gather(Fragment(r_h2o_1, m_h2o_1, isfrozen=True),
                   Fragment(r_h2o_2, m_h2o_2, isfrozen=True),
                   Fragment(None, m_mol))

    job_fde = adf_fragmentsjob(templates.singlepoint. overlay(settings), frags,
                               job_name="test_fde_fragments")

    # Perform FDE job and get dipole
    # This gets the dipole moment of the active subsystem only
    dipole_fde = run(job_fde.dipole)

    return dipole_fde
