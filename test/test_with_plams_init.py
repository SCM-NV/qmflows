import plams

plams.init()

exec(open('test/tests_using_plams/freq_ethene_DFTB.py').read())
exec(open('test/tests_using_plams/ADFGeometry_Constraints.py').read())
exec(open('test/tests_using_plams/mo_cp2k.py').read())
exec(open('test/tests_using_plams/opt_orca.py').read())
exec(open('test/tests_using_plams/linear_TS.py').read())
exec(open('test/tests_using_plams/get_properties.py').read())
exec(open('test/tests_using_plams/opt_gamess.py').read())
plams.finish()
