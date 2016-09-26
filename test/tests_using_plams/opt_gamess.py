
from nose.plugins.attrib import attr
from plams import Molecule
from qmworks import (concat, Settings, templates, zipWith)

# User Defined imports
from qmworks.packages.gamess import gamess
from qmworks.packages import run

import operator
import plams


@attr('slow')
def test_opt_gamess():
    """
    Test Optimization in Gamess using methanol in water.
    """
    methanol = Molecule('test/test_files/ion_methanol.xyz')
    methanol.properties['symmetry'] = 'Cs'

    s = Settings()
    s.specific.gamess.contrl.nzvar = 12
    s.specific.gamess.system.timlim = 2
    s.specific.gamess.system.mwords = 2
    s.specific.gamess.pcm.solvnt = 'water'
    s.specific.gamess.basis.gbasis = 'sto'
    s.specific.gamess.basis.ngauss = 3
    s.specific.gamess.guess.guess = 'huckel'
    s.specific.gamess.stapt.optol = '1d-5'
    s.specific.gamess.zmat["izmat(1)"] = "1,1,2,  1,2,3,  1,3,4,  1,3,5,  1,3,6, \n"\
                                         "2,1,2,3,  2,2,3,4,  2,2,3,5,  2,2,3,6, \n"\
                                         "3,1,2,3,4,  3,1,2,3,5,  3,1,2,3,6"

    inp = templates.geometry.overlay(s)
    methanol_geometry = gamess(inp, methanol)

    mol_opt = run(methanol_geometry.molecule)

    coords = concat([a.coords for a in mol_opt.atoms])

    expected_coords = [-0.9956983464, 0.9204754677, -0.0002616586, -0.72585581,
                       -0.0802380791, 2.18166e-05, 0.741292161, 0.0371204735,
                       -1.69738e-05, 1.1448441964, 0.5632291664, -0.9026112278,
                       1.1448447102, 0.562978981, 0.9027182521,
                       1.1454516521, -0.9993402516, 1.04943e-05]

    assert abs(sum(zipWith(operator.sub)(coords)(expected_coords))) < 1e-7

if __name__ == "__main__":
    plams.init()
    test_opt_gamess()
    plams.finish()
