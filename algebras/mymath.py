"""This module provides some basic functions and types."""
import operator
import functools
import itertools
from typing import Tuple
# TODO: use zip_longest
# TODO: reduce the use of tuple() when possible (use list comprehensions to improve the performance)


class Deg(tuple):
    """A subclass of tuple with element-wise addition and broadcast multiplication.

    All Deg instances should have the same length when added together."""
    def __new__(cls, iterable) -> "Deg":
        # noinspection PyTypeChecker
        return tuple.__new__(cls, iterable)

    def __add__(self, other):
        """Element-wise addition."""
        return Deg(map(operator.add, self, other))

    def __radd__(self, other):
        """This is implemented for supporting sum()."""
        return Deg(map(operator.add, self, other)) if other is not 0 else self

    def __sub__(self, other):
        """Element-wise addition."""
        return Deg(map(operator.sub, self, other))

    def __mul__(self, other: int) -> "Deg":
        """Broadcast multiplication."""
        return Deg(map(operator.mul, self, itertools.repeat(other)))

    def __rmul__(self, other: int):
        """Broadcast multiplication."""
        return Deg(map(operator.mul, self, itertools.repeat(other)))

    def __floordiv__(self, other) -> int:
        return min(itertools.starmap(operator.floordiv, filter(operator.itemgetter(1), zip(self, other))))


class FrozenDict(dict):
    """A subclass of dict which is hashable."""
    __setitem__ = None
    update = None
    __slots__ = "_hash",

    def __hash__(self):
        if not hasattr(self, "_hash"):
            self._hash = functools.reduce(operator.xor, map(hash, self.items()), 0)
        return self._hash


# tuple operations as monomials, everything nonnegative.
def le_tuple(t1, t2):
    """Return if t1 <= t2 element-wise."""
    return len(t1) <= len(t2) and all(map(operator.le, t1, t2))


def rstrip_tuple(t: tuple):
    """Return `t` with trailing zeroes removed."""
    if not t or t[-1]:
        return t
    right = len(t) - 1
    while right > 0 and t[right - 1] == 0:
        right -= 1
    return t[:right]


def sub_tuple(t1, t2):
    """Require le_tuple(t2, t1). Return t1 - t2 element-wise"""
    result = tuple(itertools.chain(map(operator.sub, t1, t2), t1[len(t2):]))
    return rstrip_tuple(result)


def add_tuple(t1, t2):
    """Return t1 + t2 element-wise"""
    if len(t1) < len(t2):
        return tuple(itertools.chain(map(operator.add, t1, t2), t2[len(t1):]))
    else:
        return tuple(itertools.chain(map(operator.add, t1, t2), t1[len(t2):]))


def min_tuple(t1, t2):
    return rstrip_tuple(tuple(map(min, t1, t2)))


def max_tuple(t1, t2):
    if len(t1) < len(t2):
        return tuple(itertools.chain(map(max, t1, t2), t2[len(t1):]))
    else:
        return tuple(itertools.chain(map(max, t1, t2), t1[len(t2):]))


def div_tuple(t1, t2) -> int:
    """Require le_tuple(t2, t1). Return t1 // t2."""
    return min(itertools.starmap(operator.floordiv, filter(operator.itemgetter(1), zip(t1, t2))))


def div_mod_tuple(t1, t2) -> Tuple[int, tuple]:
    """Require le_tuple(t2, t1). Return div_mod(t1, t2)."""
    q = min(itertools.starmap(operator.floordiv, filter(operator.itemgetter(1), zip(t1, t2))))
    r = tuple(itertools.chain(map(operator.sub, t1, map(operator.mul, t2, itertools.repeat(q))), t1[len(t2):]))
    return q, rstrip_tuple(r)


def add_dict(d1: tuple, d2: tuple):
    """Add tuples from tuple(sorted(d.items()))."""
    result = dict(d1)
    for gen, exp in d2:
        if gen in result:
            result[gen] += exp
        else:
            result[gen] = exp
    result = tuple(sorted(result.items()))
    return result


# binomial coefficients
def choose_mod2(m: int, n: int) -> bool:
    """Compute $\\binom{m}{n}\\text{ mod } 2$."""
    return binom_mod2(m-n, n)


def binom_mod2(m: int, n: int) -> bool:
    """Compute the binomial $(m, n)\\text{ mod } 2$."""
    return not m & n if m >= 0 and n >= 0 else 0


def multinom_mod2(*args: int) -> bool:
    """Compute the multinomial $(arg1, arg2, ...)\\text{ mod } 2$."""
    for i in args:
        if i < 0:
            return False
    s = sum(args)
    num_s = bin(s).count('1')
    num_sum = sum(bin(i).count('1') for i in args)
    return num_s == num_sum


def choose(m: int, n: int) -> int:
    """Compute $\\binom{m}{n}$."""
    n = min(n, m-n)
    if n == 0:
        return 1
    return functools.reduce(operator.mul, range(m, m - n, -1)) // functools.reduce(operator.mul, range(1, n + 1))


def binom(m: int, n: int) -> int:
    """Compute the binomial $(m,n)$."""
    return choose(m + n, n) if n >= 0 and m >= 0 else 0


# others
def xgcd(a, b):
    """Return gcd(a, b), x, y such that ax+by=gcd(a, b)."""
    x, x1 = 1, 0
    y, y1 = 0, 1
    while b:
        q, r = divmod(a, b)
        x1, x = x - q * x1, x1
        y1, y = y - q * y1, y1
        a, b = b, r
    return a, x, y


def inv_mod(a, m):
    """Return x such that ax=1 mod m."""
    d, x, y = xgcd(a, m)
    if d == 1:
        return x % m
    else:
        raise ValueError("Modular inverse does not exist.")


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


def cartanindices(k: int, n: int):
    """Return an iterator of $(i_j\\ge 0)$ such that $i_1+\\cdots+i_k=n$."""
    if k == 0:
        return
    elif k == 1:
        yield n,
        return
    for i in range(n + 1):
        for t in cartanindices(k - 1, n - i):
            yield (i,) + t


def orderedpartition(k: int, n: int):
    """Return an iterator of $(i_j>0)$ such that $i_1+\\cdots+i_k=n$."""
    if k == 0 or n < k:
        return
    elif k == 1:
        yield n,
        return
    for i in range(1, n + 1):
        for t in orderedpartition(k - 1, n - i):
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
def tex_pow(base, exp: int) -> str:
    """Return base^exp in latex."""
    if type(base) != str:
        base = str(base)
    if exp == 1:
        return base
    else:
        if "^" in base:
            base = "(" + base + ")"
        return f"{base}^{exp}" if len(str(exp)) == 1 else f"{base}^{{{exp}}}"


def tex_sub(obj, subscript) -> str:
    """Return obj_subscript in latex."""
    return f"{obj}_{subscript}" if len(str(subscript)) == 1 else f"{obj}_{{{subscript}}}"


# 73, 87, 177
