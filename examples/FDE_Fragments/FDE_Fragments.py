# Default imports

from plams import Molecule
from qmworks import (Settings, run, templates)
from qmworks.packages.SCM import adf
from qmworks.components import Fragment, adf_fragmentsjob


# ----------------------------------------------------------------

# Read the Molecule from file
m_h2o_1 = Molecule('FDE-H2O-1.xyz')
m_h2o_2 = Molecule('FDE-H2O-2.xyz')
m_mol   = Molecule('FDE-mol.xyz')

m_tot = m_mol + m_h2o_1 + m_h2o_2

settings = Settings()
settings.basis = 'SZ'
settings.specific.adf.nosymfit = ''

# Prepare first water fragment
r_h2o_1 = adf(templates.singlepoint.overlay(settings), m_h2o_1, job_name="h2o_1")

# Prepare second water fragment
r_h2o_2 = adf(templates.singlepoint.overlay(settings), m_h2o_2, job_name="h2o_2")

frozen_frags = [Fragment(r_h2o_1, [m_h2o_1]),
                Fragment(r_h2o_2, [m_h2o_2])]

job_fde = adf_fragmentsjob(templates.singlepoint. overlay(settings), m_mol, *frozen_frags)

# Perform FDE job and get dipole
# This gets the dipole moment of the active subsystem only
dipole_fde = run(job_fde.dipole)

print('FDE dipole:', dipole_fde)
