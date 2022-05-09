import numpy as np
from assertionlib import assertion

from qmflows.common import InfoMO


class TestInfoMO:
    def test_array(self) -> None:
        info_mo = InfoMO(
            eigenvalues=np.arange(3),
            eigenvectors=np.arange(12).reshape(4, 3),
        )
        ref = np.array([
            (0, [0, 3, 6, 9]),
            (1, [1, 4, 7, 10]),
            (2, [2, 5, 8, 11])
        ], dtype=[("eigenvalues", "f8"), ("eigenvectors", "f8", (4,))])

        ar = np.array(info_mo)
        assertion.eq(ar.dtype, ref.dtype)
        np.testing.assert_array_equal(ar, ref)
