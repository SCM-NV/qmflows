# from noodles import (gather, run_single)
# from pytest_mock import mocker
# from qmflows.packages.SCM import ADF_Result
# import numpy as np
# from qmflows.components.operations import (
#     find_first_job, select_max, select_min, select_max_dict, select_min_dict)


# def test_mock(mocker):
#     """
#     Mock an ADF call.
#     """
#     n = 9
#     results = [mocker.patch('qmflows.packages.SCM.ADF_Result') for _ in range(n + 1)]
#     results[n].prop = 1e3

#     rs = np.random.uniform(high=5.0, size=n)

#     for (i,), x in np.ndenumerate(rs):
#         setattr(results[i], 'prop', x)

#     wf = select_max(gather(*results), prop='prop')

#     xs = run_single(wf)

#     print(xs)

#     assert False
