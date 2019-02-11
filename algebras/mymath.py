"""This module provides some basic functions and types."""
import operator
import functools
import itertools


class Deg(tuple):
    """A subclass of tuple with element-wise addition and broadcast multiplication.

    All Deg instances should have the same length when added together."""
    def __new__(cls, iterable):
        # noinspection PyTypeChecker
        return tuple.__new__(cls, iterable)

    def __add__(self, other):
        """Element-wise addition."""
        return Deg(map(operator.add, self, other))

    def __radd__(self, other):
        """This is implemented for supporting sum()."""
        return Deg(map(operator.add, self, other)) if other is not 0 else self

    def __mul__(self, other: int):
        """Broadcast multiplication."""
        return Deg(map(lambda x: x * other, self))


class FrozenDict(dict):
    """A subclass of dict which is hashable."""
    __setitem__ = None
    update = None
    __slots__ = "_hash",

    def __hash__(self):
        if not hasattr(self, "_hash"):
            self._hash = functools.reduce(operator.xor, map(hash, self.items()), 0)
        return self._hash


# tuple operations
def leq_tuple(t1, t2):
    """Return if t1 <= t2 element-wise."""
    return len(t1) <= len(t2) and all(map(operator.le, t1, t2))


def sub_tuple(t1, t2):
    """Assert leq_tuple(t2, t1) and return t1 - t2 element-wise"""
    return tuple(itertools.chain(map(operator.sub, t1, t2), t1[len(t2):]))


def add_tuple(t1, t2):
    """Return t1 + t2 element-wise"""
    if len(t1) < len(t2):
        return tuple(itertools.chain(map(operator.add, t1, t2), t2[len(t1):]))
    else:
        return tuple(itertools.chain(map(operator.add, t1, t2), t1[len(t2):]))


def div_mod_tuple(t1, t2):
    """Assert leq_tuple(t2, t1) and return div_mod(t1, t2)."""
    q = min(map(operator.floordiv, t1, t2))
    r = tuple(itertools.chain(map(operator.sub, t1, map(lambda x: x * q, t2)), t1[len(t2):]))
    return q, r


# binomial coefficients
def choose_mod2(m: int, n: int) -> bool:
    """Compute m choose n modulo 2."""
    return binom_mod2(m-n, n)


def binom_mod2(m: int, n: int) -> bool:
    """Compute the binomial (m, n)."""
    return not m & n


def multinom_mod2(*args: int) -> bool:
    """Compute the multinomial (arg1, arg2, ...)"""
    for i in args:
        if i < 0:
            return False
    s = sum(args)
    num_s = bin(s).count('1')
    num_sum = sum(bin(i).count('1') for i in args)
    return num_s == num_sum


# others
def two_expansion(n: int):
    """If n = 2^k1 + ... + 2^kn,
    return an generator of k1, ..., kn.
    """
    k = 0
    while n:
        if n & 1:
            yield k
        n >>= 1
        k += 1


def cartan_indexing(length: int, deg: int):
    """Return an iterator of tuples for the iterated Cartan formula."""
    if length == 0:
        return
    elif length == 1:
        yield deg,
        return
    for i in range(deg + 1):
        for t in cartan_indexing(length - 1, deg - i):
            yield (i,) + t


def unique_min(iterable, *, default=None, key=None):
    """Return min if unique otherwise return None."""
    it = iter(iterable)
    minimum = next(it) if default is None else default
    unique = True
    if key is None:
        for x in it:
            if x < minimum:
                minimum = x
                unique = True
            elif x == minimum:
                unique = False
    else:
        key_minimum = key(minimum)
        for x in it:
            key_x = key(x)
            if key_x < key_minimum:
                minimum = x
                key_minimum = key_x
                unique = True
            elif key_x == key_minimum:
                unique = False
    return minimum if unique else None


# ---- latex --------
def tex_index(obj) -> str:
    """Return a string obj used to express x^obj in latex."""
    result = str(obj)
    return result if len(result) == 1 else f"{{{result}}}"


def tex_exponent(e: int) -> str:
    """Return a string ^e used in latex."""
    if e == 1:
        return ""
    else:
        return f"^{e}" if len(str(e)) == 1 else f"^{{{e}}}"


def print_tex_iter(iterable):
    for obj in iterable:
        print(f"${obj}$\\\\")

# 73, 87
