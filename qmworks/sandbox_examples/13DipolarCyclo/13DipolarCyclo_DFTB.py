# Default imports
from qmworks import Settings, templates, run, rdkitTools
from noodles import gather

# User Defined imports
from math import sqrt
from qmworks.packages.SCM import dftb
from qmworks.components import PES_scan, select_max

import plams
# ========== =============

plams.init()

hartree2kcal = 627.5095


def bond_distance(r1, r2):
    return sqrt(sum((x - y) ** 2 for x, y in zip(r1, r2)))

startvalue = {'C': 2.1, 'N': 1.9, 'O': 1.8}

settings = Settings()
settings.specific.dftb.dftb.scc

ethene = rdkitTools.smiles2plams('C=C')
E_ethene = run(dftb(templates.geometry.overlay(settings), ethene)).energy

molecules = set()
result = {}
for a1 in ('C', 'N', 'O'):
    for a2 in ('N', 'O', 'S'):
        for a3 in ('C', 'N', 'O'):
            # reactant
            smiles = 'C1' + a1 + a2 + a3 + 'C1'
            reverse = 'C1' + a3 + a2 + a1 + 'C1'
            if reverse in molecules:  # skip duplicates
                continue
            molecules.add(smiles)
            print(smiles)
            product  = rdkitTools.smiles2plams(smiles)
            E_product = dftb(templates.geometry.overlay(settings),
                             product).energy

            smiles = {'C': '[CH2]=', 'N': '[NH]=', 'O': 'O='}[a1] +\
                     {'N': '[NH+]', 'O': '[O+]', 'S': '[S+]'}[a2] +\
                     {'C': '[CH2-]', 'N': '[NH-]', 'O': '[O-]'}[a3]
            reactant = rdkitTools.smiles2plams(smiles)
            E_reactant = dftb(templates.geometry.overlay(settings),
                              reactant).energy

            constraint1 = "dist 1 2"
            constraint2 = "dist 4 5"

            sv1 = startvalue[a1]
            sv2 = startvalue[a3]
            # scan input
            if a1 == a3:
                # symmetric molecule
                scan = {'constraint': [constraint1, constraint2],
                        'surface': {'nsteps': 4, 'start': [sv1, sv2],
                                    'stepsize': [0.1, 0.1]}}
            else:
                scan = {'constraint': constraint1,
                        'surface': {'nsteps': 4, 'start': sv1, 'stepsize': 0.1},
                        'scan': {'constraint': constraint2,
                                 'surface': {'nsteps': 4, 'start': sv2,
                                             'stepsize': 0.1}
                                 }
                        }

            # returns a set of results object containing the output of
            # each point in the scan
            LT = PES_scan(dftb, settings, product, scan)
                        
            # Gets the object presenting the molecule
            # with the maximum energy calculated from the scan
            apprTS = select_max(LT, 'energy')
            
            # Run the TS optimization, using the default TS template
            TS = dftb(templates.ts.overlay(settings), apprTS.molecule)

            # Actual execution of the jobs
            reactives = [E_reactant, E_product, TS]
            E_reactant, E_product, TS = run(gather(*reactives), n_processes=1)

            # Retrieve the molecular coordinates
            mol = TS.molecule
            d1 = bond_distance(mol.atoms[0].coords, mol.atoms[1].coords)
            d2 = bond_distance(mol.atoms[3].coords, mol.atoms[4].coords)
            
            Eact = (TS.energy - E_ethene - E_reactant) * hartree2kcal
            Ereact = (E_product - E_ethene - E_reactant) * hartree2kcal
            result[smiles] = [Eact, Ereact, d1, d2]
 
            print("Reactant                Eact  Ereact   Bond1   Bond2")
            for smiles in result:
                print("{0:20s} {1:7.1f} {2:7.1f} {3:7.2f} \
                {4:7.2f}".format(smiles, *result[smiles]))
