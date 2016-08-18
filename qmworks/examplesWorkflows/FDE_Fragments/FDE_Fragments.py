# Default imports
from qmworks import Settings, templates, run
from noodles import gather, schedule, Storable
from plams import Molecule

# User Defined imports
from qmworks.packages.SCM import adf
# from qmworks.components import adffragmentsjob

import plams
# ========== =============


class Fragment(Storable):
    def __init__(self, result, mol_list):
        self.result = result
        self.mol_list = gather(*mol_list)

# This function should be made more generic and go to "components"
# ----------------------------------------------------------------


@schedule
def fragmentsjob(settings, mol, *frozen_frags):
    mol_tot = Molecule()
    frag_settings = Settings()
    for i, frag in enumerate(frozen_frags):
        i += 1
        frag_id = 'frag' + str(i)
        for m in frag.mol_list:
            for a in m:
                a.fragment = frag_id
            mol_tot += m
        path = frag.result.archive.path + ' type=FDE'
        frag_settings.specific.adf.fragments[frag_id] = path
        frag_settings.specific.adf.fde.PW91k = ""
    mol_tot += mol

    return adf(settings.overlay(frag_settings), mol_tot)
# ----------------------------------------------------------------


plams.init()

# Read the Molecule from file
m_h2o_1 = Molecule('FDE-H2O-1.xyz')
m_h2o_2 = Molecule('FDE-H2O-2.xyz')
m_mol   = Molecule('FDE-mol.xyz')

m_tot = m_mol + m_h2o_1 + m_h2o_2

print(m_tot)

settings = Settings()
settings.basis = 'SZ'
settings.specific.adf.nosymfit = ''

# Prepare first water molecule
r_h2o_1 = adf(templates.singlepoint.overlay(settings), m_h2o_1)

# Prepare second water molecule
r_h2o_2 = adf(templates.singlepoint.overlay(settings), m_h2o_2)

frozen_frags = [Fragment(r_h2o_1, [m_h2o_1]),
                Fragment(r_h2o_2, [m_h2o_2])]

fde_res = run(fragmentsjob(settings, m_mol, *frozen_frags))

print(fde_res.dipole)


