"""
Provides the basic functions
"""


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
    """ return a string that is to used to express a^n in latex """
    result = str(obj)
    return result if len(result) == 1 else "{" + result + "}"

# 54
