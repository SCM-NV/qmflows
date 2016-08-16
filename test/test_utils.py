
from qmworks.utils import (concat, concatMap)


def test_concat():
    """
    Test list concatenation
    """
    xss = [[1], [2]]
    assert concat(xss) == [1, 2]


def test_concatMap():
    """
    concatMap f == concat . map f
    """
    f = lambda x: [x * 2]
    xs = [1, 2]
    assert concatMap(f, xs) == list(concat(map(f, xs)))
