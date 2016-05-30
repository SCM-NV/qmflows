# from qmworks import Settings, templates
# from qmworks.packages import registry
# import plams
# from noodles import run_fireworks, gather

# # User Defined imports
# from qmworks.packages.SCM import adf


# def test_ADFGeometry_Constraint():
#     """
#     Test "freeze" and "selected_atoms" keywords for constrained geometry optimizations
#     """
#     plams.init()
#     an = plams.Molecule('test_files/an.xyz', 'xyz')
    
#     # optimize only H atoms
#     s = Settings()
#     s.freeze = [1, 2, 3]
#     result1 = adf(templates.geometry.overlay(s), an)
#     geom1 = result1.molecule

#     r = run_fireworks(geom1, remote_db="felipe@145.100.59.99", registry=registry)
#     print(r)
#     assert False
#     # # optimize only H atoms
#     # s = Settings()
#     # s.selected_atoms = ['H']
#     # result2 = adf(templates.geometry.overlay(s), an)
#     # geom2 = result2.molecule
    
#     # r = run_process(gather(geom1, geom2), n_processes = 1, registry=registry, init=init, finish=finish)
#     # #print(str(r[0]), str(r[1]))
#     # assert str(r[0]) == str(r[1])
