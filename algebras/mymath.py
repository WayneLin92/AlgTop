""" Provides some basic functions and types """
import operator


class Deg(tuple):
    """ A simple version of numpy.array for modeling multi-degrees """
    def __new__(cls, iterable):
        # noinspection PyTypeChecker
        return tuple.__new__(cls, iterable)

    def __add__(self, other):
        """ element-wise addition """
        return Deg(map(operator.add, self, other))

    def __radd__(self, other):
        """ this is implemented for supporting sum() """
        return Deg(map(operator.add, self, other)) if other is not 0 else self

    def __mul__(self, other):
        """ assert type(other) is int """
        return Deg(map(lambda x: x * other, self))


def choose_mod2(m: int, n: int) -> bool:
    """ Compute m choose n modulo 2~ """
    return binom_mod2(m-n, n)


def binom_mod2(m: int, n: int) -> bool:
    return not m & n


def multi_nom_mod2(*args: int) -> bool:
    for i in args:
        if i < 0:
            return 0
    s = sum(args)
    num_s = bin(s).count('1')
    num_sum = sum(bin(i).count('1') for i in args)
    return num_s == num_sum


def two_expansion(n: int):
    """
    If n = 2^k1 + ... + 2^kn,
    return an iterator of k1, ..., kn.
    """
    k = 0
    while n:
        if n & 1:
            yield k
        n >>= 1
        k += 1


def cartan_indexing(length, deg):
    """ return an iterator of tuples for the Cartan formula """
    if length == 0:
        return
    elif length == 1:
        yield deg,
        return
    for i in range(deg + 1):
        for t in cartan_indexing(length - 1, deg - i):
            yield (i,) + t


def tex_index(obj) -> str:
    """ return a string that is to used to express x^obj in latex """
    result = str(obj)
    return result if len(result) == 1 else "{" + result + "}"

# 73
