"""Test utility functions."""
from pymonad import curry
from qmflows.utils import zipWith


def test_zipwith():
    """zipWith f xs ys ==  map (f) (zip  xs ys)."""
    @curry
    def f(x, y):
        return x * y

    xs = range(3)
    ys = range(3, 6)

    assert zipWith(f)(xs)(ys) == [f(x, y) for (x, y) in zip(xs, ys)]
