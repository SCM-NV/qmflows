from pymonad import curry
from qmflows.utils import zipWith, zipWith3


def test_zipwith():
    """zipWith f xs ys ==  map (f) (zip  xs ys)."""
    @curry
    def f(x, y):
        return x * y

    xs = range(3)
    ys = range(3, 6)

    assert zipWith(f)(xs)(ys) == [f(x, y) for (x, y) in zip(xs, ys)]


def test_zipwith3():
    """zipWith f xs ys zs ==  map (f) (zip  xs ys zs)."""
    @curry
    def f(x, y, z):
        return x * y + z

    xs = range(3)
    ys = range(3, 6)
    zs = range(6, 9)

    assert zipWith3(f)(xs)(ys)(zs) == [f(x, y, z)
                                       for (x, y, z) in zip(xs, ys, zs)]
