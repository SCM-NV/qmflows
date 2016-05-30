
__author__ = "Felipe Zapata"

__all__ = ['floatArray', 'applyN', 'calculateUniqueLabel', 'chunksOf',
           'concat', 'concatMap', 'fst', 'flatten', 'head', 'headTail', 'product',
           'repeatN', 'replicate', 'snd', 'str2Float', 'str2Int', 'swapCoeff', 'zipWith',
           'zipWith3', 'zipWith4', 'odd', 'even', 'fac', 'binomial', 'settings2Dict',
           'dict2Setting', 'stringVal2Dict', 'stringDict2Settings', 'recursiveLookup']

# ================> Python Standard  and third-party <==========
import numpy    as np
import operator as op

from functools import reduce
from itertools import chain
from pymonad   import curry

# ======================> List Functions <========================


def applyN(f, n, x):
    """
    apply a function n times, returning the results as a list
    """
    acc = []
    for i in range(n):
        r = f(x)
        acc.append(r)
    return acc


def calculateUniqueLabel(xs):
    """Return a list of unique tokens"""
    rs = set()
    for x in xs:
        if x not in rs:
            rs.add(x)
    return rs


def chunksOf(xs, n):
    """Yield successive n-sized chunks from xs"""
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def concat(xss):
    """The concatenation of all the elements of a list"""
    return list(chain(*xss))


def concatMap(f, xss):
    """Map a function over all the elements of a container and concatenate the resulting lists"""
    return concat(list(map(f, xss)))


def fst(xs):
    return xs[0]


def flatten(xs):
    return reduce(lambda x, y: x + y, xs)


def head(xs):
    return xs[0]


def headTail(xs):
    it = iter(xs)
    head = next(it)
    tail = list(it)
    return (head, tail)


def product(xs):
    return reduce(op.mul, xs)


def repeatN(n, a):
    for x in range(n):
        yield a


def replicate(n, a):
    return list(repeatN(n, a))


def snd(xs):
    return xs[1]


def splitAtNone(xs):
    """Yield Elements of the list until one of the elements is None"""
    for i, x in enumerate(xs):
        if x is None:
            break
    return xs[:i]


def str2Float(xs):
    return [float(x) for x in xs.split()]


def str2Int(xs):
    return [int(x) for x in xs.split()]


@curry
def swapCoeff(n, rss):

    if n == 1:
        return [rss]
    else:
        return [[rs[i::n] for i in range(n)] for rs in rss]

    
@curry
def zipWith(f, xs, ys):
    """zipWith generalises zip by zipping with the function given as the first argument"""
    return [f(*rs) for rs in zip(xs, ys)]


@curry
def zipWith3(f, xs, ys, zs):
    """
    The zipWith3 function takes a function which combines three elements,
    as well as three lists and returns a list of their point-wise combination.
    """
    return [f(*rs) for rs in zip(xs, ys, zs)]


@curry
def zipWith4(f, ws, xs, ys, zs):
    return [f(*rs) for rs in zip(ws, xs, ys, zs)]


# ============> Integer functions <================
def odd(x):
    return x % 2 != 0


def even(x):
    return not(odd(x))


def fac(x):
    if x == 0:
        return 1
    else:
        return float(product(range(1, x + 1)))


def binomial(n, k):
    if k == n:
        return 1.0
    elif k >= 0 and k < n:
        return fac(n) / (fac(k) * fac(n - k))
    else:
        return 0.0


# ================> Dict Functions
from qmworks.settings   import Settings


def settings2Dict(s):
    """
    Transform a Settings object into a dict
    """
    d = {}
    for k, v in s.items():
        if not isinstance(v, Settings):
            d[k] = v
        else:
            d[k] = settings2Dict(v)

    return d


def dict2Setting(d):
    """
    """
    r = Settings()
    for k, v in d.items():
        if isinstance(v, dict):
            r[k] = dict2Setting(v)
        else:
            r[k] = v

    return r


def stringVal2Dict(xs, v):
    """
    """
    ys = reversed(xs.split('.'))

    go = lambda acc, x: {x: acc}
    return  reduce(go, ys, v)


def stringDict2Settings(xs):
    r = Settings()
    for k, v in xs.items():
        d = stringVal2Dict(k, v)
        s = dict2Setting(d)
        r = r.merge(s)
    return r


def recursiveLookup(d, k):
    v = None
    try:
        v = d[k]
    except KeyError:
        for x in d.values():
            if isinstance(x, dict):
                v = recursiveLookup(x, k)

    return v


# -------------Array Functions -----------------------#
floatArray = np.vectorize(float)


