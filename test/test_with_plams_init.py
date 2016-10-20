from nose.plugins.attrib import attr
import importlib
import plams


# @attr('slow')
# def test_packages():
#     """
#     ``modules`` is a dictionary contains as key the python modules names and
#     a list of functions as values.
#     """
#     modules = {
#         'freq_ethene_DFTB': ['test_freq'],
#         'ADFGeometry_Constraints': ['test_ADFGeometry_Constraint'],
#         'mo_cp2k': ['test_ethylene'],
#         'opt_orca': ['test_opt_orca', 'test_methanol_opt_orca'],
#         'linear_TS': ['test_linear_ts'],
#         'get_properties': ['test_dftb_props', 'test_adf_props'],
#         'opt_gamess': ['test_opt_gamess']
#     }

#     plams.init()

#     root_folder = 'test.tests_using_plams'

#     # Load all the tests modules
#     for module_name, fs in modules.items():
#         print("testing modules: ", module_name)
#         mod = root_folder + '.' + module_name
#         m = importlib.import_module(mod)
#         # Execute each test function inside the module
#         for fun in fs:
#             getattr(m, fun)()

#     plams.finish()
