"""This module provides some basic functions and types."""
import operator
from functools import reduce
from itertools import repeat, starmap, chain, zip_longest
from typing import Tuple, Iterable


class Vector(tuple):
    """A subclass of tuple with element-wise addition and scalar multiplication.

    All Vector instances should have the same length when added together."""
    __slots__ = []

    def __new__(cls, iterable) -> "Vector":
        # noinspection PyTypeChecker
        return tuple.__new__(cls, iterable)

    def __add__(self, other):
        """Element-wise addition."""
        return Vector(map(operator.add, self, other))

    def __iadd__(self, other):
        """Element-wise addition."""
        return Vector(map(operator.add, self, other))

    def __radd__(self, other):
        """This is implemented for supporting sum()."""
        return Vector(map(operator.add, self, other)) if other != 0 else self

    def __sub__(self, other):
        """Element-wise addition."""
        return Vector(map(operator.sub, self, other))

    def __isub__(self, other):
        """Element-wise addition."""
        return Vector(map(operator.sub, self, other))

    def __mul__(self, other) -> "Vector":
        """Scalar multiplication."""
        return Vector(map(operator.mul, self, repeat(other)))

    def __rmul__(self, other):
        """Scalar multiplication."""
        return Vector(map(operator.mul, self, repeat(other)))

    def __floordiv__(self, other):
        """Floor division by a Vector or a scalar."""
        if type(other) is Vector:
            return min(starmap(operator.floordiv, filter(operator.itemgetter(1), zip(self, other))))
        else:
            return Vector(map(operator.floordiv, self, repeat(other)))

    def __truediv__(self, other) -> "Vector":
        """Division by a scalar."""
        return Vector(map(operator.truediv, self, repeat(other)))


# tuple operations especially for dense monomials
def le_tuple(t1, t2):
    """Return if t1_i <= t2_i."""
    return len(t1) <= len(t2) and all(map(operator.le, t1, t2))


def rstrip_tuple(t: tuple):
    """Remove trailing zeroes in `t`."""
    if not t or t[-1]:
        return t
    right = len(t) - 1
    while right > 0 and t[right - 1] == 0:
        right -= 1
    return t[:right]


def sub_tuple(t1, t2):
    """Require len(t2)<len(t1). Return t1 - t2 element-wise."""
    result = tuple(chain(map(operator.sub, t1, t2), t1[len(t2):]))
    return rstrip_tuple(result)


def add_tuple(t1, t2):
    """Return t1 + t2 element-wise"""
    return tuple(starmap(operator.add, zip_longest(t1, t2, fillvalue=0)))


def mul_tuple(t, k):
    """Return Scalar product t * k."""
    return tuple(map(operator.mul, t, repeat(k)))


def min_tuple(t1, t2):
    """return (min(t1_i, t2_i), ...)."""
    return rstrip_tuple(tuple(map(min, t1, t2)))


def max_tuple(t1, t2):
    """return (max(t1_i, t2_i), ...)."""
    return tuple(starmap(max, zip_longest(t1, t2, fillvalue=0)))


def div_tuple(t1, t2) -> int:
    """Require le_tuple(t2, t1). Return the largest q such that t2 * q <= t1."""
    return min(starmap(operator.floordiv, filter(operator.itemgetter(1), zip(t1, t2))))


def div_mod_tuple(t1, t2) -> Tuple[int, tuple]:
    """Require le_tuple(t2, t1). Return div_mod(t1, t2)."""
    q = min(starmap(operator.floordiv, filter(operator.itemgetter(1), zip(t1, t2))))
    r = tuple(chain(map(operator.sub, t1, map(operator.mul, t2, repeat(q))), t1[len(t2):]))
    return q, rstrip_tuple(r)


# tuple operations especially for sparse monomials
def add_dtuple(d1, d2):
    """Return d1 + d2 as sparse vectors."""
    result = dict(d1)
    for gen, exp in d2:
        if gen in result:
            result[gen] += exp
        else:
            result[gen] = exp
    return tuple(sorted(result.items()))


def sub_dtuple(d1, d2):
    """Return d1 - d2 as sparse vectors."""
    result = dict(d1)
    for gen, exp in d2:
        if gen in result:
            result[gen] -= exp
        else:
            result[gen] = -exp
    return tuple(sorted((k, v) for k, v in result.items() if v))


def le_dtuple(d1, d2):
    """Return if d1_i <= d2_i as sparse vectors."""
    d2_dict = dict(d2)
    return all(gen in d2_dict and exp <= d2_dict[gen] for gen, exp in d1)


def min_dtuple(d1, d2):
    """return (min(d1_i, d2_i), ...)."""
    d1_dict = dict(d1)
    result = {}
    for gen, exp in d2:
        if gen in d1_dict:
            result[gen] = min(exp, d1_dict[gen])
    return tuple(sorted(result.items()))


def max_dtuple(d1, d2):
    """return (max(d1_i, d2_i), ...)."""
    result = dict(d1)
    for gen, exp in d2:
        result[gen] = max(exp, result[gen]) if gen in result else exp
    return tuple(sorted(result.items()))


def div_dtuple(d1, d2) -> int:
    """Require le_dtuple(t2, t1). Return the largest q such that d2 * q <= d1"""
    d1_dict = dict(d1)
    return min(d1_dict[gen] // exp for gen, exp in d2)


def div_mod_dtuple(d1, d2) -> Tuple[int, tuple]:
    """Require le_dtuple(d2, d1). Return div_mod(d1, d2)."""
    d1_dict = dict(d1)
    q = min(d1_dict[gen] // exp for gen, exp in d2)
    for gen, exp in d2:
        d1_dict[gen] -= exp * q
    return q, tuple(sorted((k, v) for k, v in d1_dict.items() if v))


# binomial coefficients
def choose_mod2(m: int, n: int) -> bool:
    R"""Compute $\binom{m}{n}\text{ mod } 2$."""
    return binom_mod2(m-n, n)


def binom_mod2(m: int, n: int) -> bool:
    R"""Compute the binomial $(m, n)\text{ mod } 2$."""
    return m >= 0 and n >= 0 and not m & n


def multinom_mod2(*args: int) -> bool:
    R"""Compute the multinomial $(arg1, arg2, ...)\text{ mod } 2$."""
    s = 0
    return all(i >= 0 and not i & ((s := s + i) - i) for i in args)


def choose(m: int, n: int) -> int:
    R"""Compute $\binom{m}{n}$."""
    n = min(n, m-n)
    if n == 0:
        return 1
    return reduce(operator.mul, range(m, m - n, -1)) // reduce(operator.mul, range(1, n + 1))


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
    R"""Return an iterator of $(i_j\ge 0)$ such that $i_1+\cdots+i_k=n$."""
    if k == 0:
        return
    elif k == 1:
        yield n,
        return
    for i in range(n + 1):
        for t in cartanindices(k - 1, n - i):
            yield (i,) + t


def orderedpartition(k: int, n: int):
    R"""Return an iterator of $(i_j>0)$ such that $i_1+\cdots+i_k=n$."""
    if k == 0 or n < k:
        return
    elif k == 1:
        yield n,
        return
    for i in range(1, n + 1):
        for t in orderedpartition(k - 1, n - i):
            yield (i,) + t


def unique_min(iterable: Iterable, *, default=None, key=None):
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


def prod_algs(iterable: Iterable, default=None):
    """Return prod of elements in algebras."""
    it = iter(iterable)
    initial = next(it) if default is None else default
    result = initial.data
    for x in it:
        result = x.mul_data(result, x.data)
    return type(initial)(result)


def clip(x, min_, max_):
    """Clip value `x` by [min_, max_]."""
    return min_ if x < min_ else (max_ if x > max_ else x)


def interpolation(alpha, t1: tuple, t2: tuple):
    """Return alpha * t1 + (1 - alpha * t2)."""
    return tuple(alpha * i + (1 - alpha) * j for i, j in zip(t1, t2))


def get_from_singleton(singleton):
    """When singleton is a set of one element, return this element."""
    for e in singleton:
        return e


# ---- latex --------
def tex_pow(base, exp: int) -> str:
    """Return base^exp in latex."""
    if type(base) != str:
        base = str(base)
    if exp == 1:
        return base
    else:
        if tex_outside_delimiter(base, "^"):
            base = "(" + base + ")"
        return f"{base}^{exp}" if len(str(exp)) == 1 else f"{base}^{{{exp}}}"


def tex_sub(obj, subscript) -> str:
    """Return obj_subscript in latex."""
    return f"{obj}_{subscript}" if len(str(subscript)) == 1 else f"{obj}_{{{subscript}}}"


def tex_parenthesis(obj):
    """Return obj with parenthesis if there is a plus or minus sign."""
    result = str(obj)
    return f"({result})" if "+" in result or "-" in result else result


def tex_braces(obj):
    """Return obj with curly braces if there are more than one character."""
    result = str(obj)
    return f"{{{result}}}" if len(result) > 1 else result


def tex_outside_delimiter(text: str, symbol: str):
    """Return if symbol appears in text and is outside any pair of delimiters including ()[]{}"""
    left, right = "([{", ")]}"
    left_minus_right = 0
    for c in text:
        if c in left:
            left_minus_right += 1
        elif c in right:
            left_minus_right -= 1
        elif c == symbol and left_minus_right == 0:
            return True
    return False
