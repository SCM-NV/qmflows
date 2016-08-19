import plams

plams.init()

exec(open('test/tests_using_plams/freq_ethene.py').read())
exec(open('test/tests_using_plams/ADFGeometry_Constraints.py').read())
exec(open('test/tests_using_plams/mo_cp2k.py').read())
exec(open('test/tests_using_plams/opt_orca.py').read())

plams.finish()
