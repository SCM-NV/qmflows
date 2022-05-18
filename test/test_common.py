import numpy as np
from assertionlib import assertion

from qmflows.common import CP2KInfoMO


class TestInfoMO:
    def test_array(self) -> None:
        info_mo = CP2KInfoMO(
            eigenvalues=np.arange(3),
            eigenvectors=np.arange(12).reshape(4, 3),
            orb_index=np.array([5, 6, 7]),
            occupation=np.array([2.0, 2.0, 0.0]),
        )
        ref = np.array([
            (0, [0, 3, 6, 9], 5, 2),
            (1, [1, 4, 7, 10], 6, 2),
            (2, [2, 5, 8, 11], 7, 0),
        ], dtype=[
            ("eigenvalues", "f8"),
            ("eigenvectors", "f8", (4,)),
            ("orb_index", "i8"),
            ("occupation", "f8"),
        ])

        ar = np.array(info_mo)
        assertion.eq(ar.dtype, ref.dtype)
        np.testing.assert_array_equal(ar, ref)
