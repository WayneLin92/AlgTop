from typing import Set, Tuple, Optional, Union

import algebras.BaseAlgebras as BC
from algebras.mymath import binom_mod2
from algebras.operations import DyerLashof


class MyDyerLashof(BC.OperationsMod2):
    # -- Algebra -------------------------
    def __init__(self, data: Union[DyerLashof, set, tuple]):
        if type(data) is DyerLashof:
            self.data = set()
            for m in data.data:
                self.data ^= Q2sQ(m).data
        else:
            super().__init__(data)

    @staticmethod
    def str_mon(mon: tuple) -> str:
        result = ""
        for i in mon:
            if i >= 10:
                result += "\sQ^{{{}}}".format(i)
            else:
                result += "\sQ^{}".format(i)
        if result == "":
            result = "1"
        return result

    # -- OperationsMod2 ------------------
    def __mul__(self, other) -> "MyDyerLashof":
        if not isinstance(other, MyDyerLashof):
            return NotImplemented
        else:
            return BC.AlgebraMod2.__mul__(self, other)

    def simplify(self, degree: Optional[int]=None):
        data = set()
        for m in self.data:
            if self.is_admissible(m)[0]:
                data ^= {m}
            else:
                data ^= simplify_sQ(m)
        self.data = data
        return self

    # the following three methods are for testing for now
    @staticmethod
    def is_null(mon: tuple, degree=None) -> bool:
        return False

    @staticmethod
    def is_admissible(mon: tuple) -> Tuple[bool, Optional[tuple]]:
        # ordering (0)
        for i in range(len(mon) - 1):
            if mon[i] > mon[i + 1]:
                return False, (i, 0)
        # repeating (1) and admissibility (2)
        for i in range(len(mon) - 1):
            if mon[i] == mon[i + 1] > 0:
                return False, (i, 1)
            elif 2 * mon[i] > mon[i + 1]:
                return False, (i, 2)
        # divisibility (3)
        for i in range(len(mon)-1, 0, -1):
            if mon[i] % (1 << i) != 0:
                return False, (i, 3)
        return True, None

    @staticmethod
    def adem(mon: tuple, index: tuple) -> Set[tuple]:
        if index[1] == 0:
            return {tuple(sorted(mon))}
        elif index[1] == 1:
            i = index[0]
            return {tuple(sorted(mon[0:i] + mon[i+2:] + (0, mon[i] * 2)))}
        elif index[1] == 2:
            i = index[0]
            r, s = mon[i], mon[i+1]
            return {mon[0:i] + m + mon[i+2:] for m in _my_adem(r, s)}
        elif index[1] == 3:
            i = index[0]
            n = 0
            while mon[i] % (1 << n) == 0:
                n += 1
            I, s = mon[:n], mon[i]
            return {tuple(I[j] - J[j] for j in range(n)) + mon[n:i] + (s + k,) + mon[i+1:]
                    for k in range(1 << n - 1, sum(I) + 1, 1 << n - 1) if binom_mod2(k + (1 << n - 1), s)
                    for J in _indexing_set(I, k) if _is_good(I)}  # ##################

    @classmethod
    def gen(cls, *n) -> "MyDyerLashof":
        return cls({n})

    # methods ---------------------
    def convert_to_dyerlashof(self):
        return sum((sQ2Q(m) for m in self.data), DyerLashof.zero())


def _my_adem(r: int, s: int) -> Set[tuple]:
    """ The Adem relation for \sR. Assert r <= s """
    assert r <= s
    if s % 2:
        return {(r-k, s+k) for k in range(1, r+1, 2) if binom_mod2(k+1, s)}
    if r % 2:
        result = set()
        for r1, s1 in _my_adem(2 * r, 2 * s):
            result ^= {(r1 >> 1, s1 >> 1)} if not s1 & 2 else _my_adem(r1 >> 1, s1 >> 1)
        return result
    n = r.bit_length() - 1
    if n == s.bit_length() - 1:
        # s = 2^n - 2^m + ...
        m = n
        while binom_mod2(1 << m-1, s) == 0:
            m -= 1
        if r < (1 << n) + (1 << m-1):
            set1 = {(r - (1 << n) + k, s + (1 << n) - k) for k in range(2, (1 << n) - (1 << (m - 1)) + 2, 2) if
                    binom_mod2(k, r - 1)}
            set2 = {(s - (1 << n) + k, r + (1 << n) - k) for k in range(2, (1 << m-1) + 2, 2) if binom_mod2(k, s - 1)}
            return set1 ^ set2 ^ {(0, r + s)}
        else:
            set1 = {(r - (1 << n) + k, s + (1 << n) - k) for k in range(2, (1 << n + 1) - r, 2) if
                    binom_mod2(k, r - 1)}
            set2 = {(s - (1 << n) + k, r + (1 << n) - k) for k in range(2, (1 << n + 1) - s, 2) if
                    binom_mod2(k, s - 1)}
            return set1 ^ set2 ^ {(0, r + s)}
    else:
        if binom_mod2(1 << n, s) == 0:  # s = 2^{n+1} + 2^n + ...
            return {(r1 + (1 << n), s1 + (1 << n+1)) for r1, s1 in _my_adem(r - (1 << n), s - (1 << n + 1))}
        elif binom_mod2(1 << n-1, r):  # r = 2^n + r_1, r_1 < 2^{n-2}
            set1 = _my_adem(r - (1 << n - 1), s - (1 << n))  # recursive
            set2 = {(r1, s1 + (3 << n-1)) for r1, s1 in set1} ^ \
                   {(r1 + (1 << n - 1), s1 + (1 << n)) for r1, s1 in set1} ^ \
                   {(r - (1 << n - 1), s + (1 << n-1))}
            return set2
        else:  # r = 2^n + 2^{n-1} + ..., s = 2^{n+1} + s_1 + ..., s_1 < 2^n
            if s < (0b101 << n-1):
                set1 = {(r - (1 << n) + k, s + (1 << n) - k) for k in range(2, (1 << n+1) - r, 2) if
                        binom_mod2(k, r-1)}
                set2 = {(s - (1 << n + 1) + k, r + (1 << n + 1) - k) for k in range(2, (1 << n) * 3 - s, 2) if
                        binom_mod2(k, s-1)}
                return set1 ^ set2 ^ {(s - (1 << n), r + (1 << n))}
            else:  # r = 2^n + 2^{n-1} + ..., s = 2^{n+1} + 2^{n-1} + ...
                if r - (1 << n) <= s - (1 << n + 1):
                    # s = 2^{n+1} + 2^n - 2^m + ...
                    m = n - 1
                    while (1 << m - 1) & s:
                        m -= 1
                    if r < (3 << n - 1) + (1 << m - 1):
                        set1 = {(r - (1 << n) + k, s + (1 << n) - k) for k in
                                range((1 << n - 1) - (1 << m - 1) + 2, (1 << n + 1) - r + 2, 2)
                                if binom_mod2(k, r - 1)}
                        set2 = {(s - (1 << n + 1) + k, r + (1 << n + 1) - k) for k in
                                range(2, (1 << m - 1) + 2, 2) if binom_mod2(k, s - 1)}
                        set3 = {(r - (1 << n - 1) + k, s + (1 << n - 1) - k) for k in
                                range(2, (1 << n - 1) - (1 << m - 1) + 2, 2) if binom_mod2(k, r - 1)}
                        set4 = {(s - (3 << n - 1) + k, r + (3 << n - 1) - k) for k in
                                range(2, (1 << m - 1) + 2, 2) if binom_mod2(k, s - 1)}
                        return set1 ^ set2 ^ set3 ^ set4
                    else:
                        set2 = {(s - (1 << n + 1) + k, r + (1 << n + 1) - k) for k in
                                range(2, (3 << n) - s + 2, 2) if binom_mod2(k, s - 1)}
                        set3 = {(r - (1 << n - 1) + k, s + (1 << n - 1) - k) for k in
                                range(2, (1 << n + 1) - r + 2, 2) if binom_mod2(k, r - 1)}
                        set4 = {(s - (3 << n - 1) + k, r + (3 << n - 1) - k) for k in
                                range(2, (3 << n) - s + 2, 2) if binom_mod2(k, s - 1)}
                        return set2 ^ set3 ^ set4
                else:
                    # r = 2^{n+1}-2^m
                    m = n - 1
                    while (1 << m - 1) & r:
                        m -= 1
                    if s < (0b101 << n - 1) + (1 << m - 1):
                        set1 = {(r - (1 << n) + k, s + (1 << n) - k) for k in
                                range((1 << m - 1) + 2, (1 << n + 1) - r + 2, 2) if binom_mod2(k, r - 1)}
                        set2 = {(s - (1 << n + 1) + k, r + (1 << n + 1) - k) for k in
                                range(2, (1 << n - 1) - (1 << m - 1) + 2, 2) if binom_mod2(k, s - 1)}
                        set3 = {(r - (1 << n - 1) + k, s + (1 << n - 1) - k) for k in
                                range(2, (1 << m - 1) + 2, 2) if binom_mod2(k, r - 1)}
                        set4 = {(s - (3 << n - 1) + k, r + (3 << n - 1) - k) for k in
                                range(2, (1 << n - 1) - (1 << m - 1) + 2, 2) if binom_mod2(k, s - 1)}
                        return set1 ^ set2 ^ set3 ^ set4
                    else:
                        set2 = {(s - (1 << n + 1) + k, r + (1 << n + 1) - k) for k in
                                range(2, (3 << n) - s + 2, 2) if binom_mod2(k, s - 1)}
                        set3 = {(r - (1 << n - 1) + k, s + (1 << n - 1) - k) for k in
                                range(2, (1 << n + 1) - r + 2, 2) if binom_mod2(k, r - 1)}
                        set4 = {(s - (3 << n - 1) + k, r + (3 << n - 1) - k) for k in
                                range(2, (3 << n) - s + 2, 2) if binom_mod2(k, s - 1)}
                        return set2 ^ set3 ^ set4


def _indexing_set(I, k):
    """ iterator of J such that deg(J)=k and J<=I"""
    if len(I) == 1:
        if k <= I[0]:
            yield k,
        return
    for jn in range(min(I[-1], k) + 1):
        for J in _indexing_set(I[:-1], k - jn):
            yield J + (jn,)


def _is_good(I):
    """ return if I shows up in seq_minus """
    d = sum(I)
    n = len(I)
    m = 0
    while (1 << m + n - 1) < (1 << n) + 3 * d:
        m += 1
    s = (1 << m + n - 1) - (1 << n) - d
    return (0,) * n + (d + s,) in simplify_sQ(I + (s,))


# ------------------ convert sQ^I to Q^J --------------------------------
def PQ_mon(r: int, m: tuple) -> DyerLashof:
    """ return P^rQ^m x assuming P^r x=0 for r > 0 """
    if len(m) == 0:
        return DyerLashof.unit() if r == 0 else DyerLashof.zero()
    elif len(m) == 1:
        return DyerLashof.gen(m[0] - r) if binom_mod2(r, m[0] - 2 * r) else DyerLashof.zero()
    else:
        s = m[0]  # recursive
        return sum((DyerLashof.gen(s - r + i) * PQ_mon(i, m[1:])
                    for i in range(max(0, r - s // 2), r // 2 + 1)
                    if binom_mod2(r - 2 * i, s - 2 * r + 2 * i)), DyerLashof.zero())


def PQ(r: int, q: DyerLashof) -> DyerLashof:
    return sum((PQ_mon(r, i) for i in q.data), DyerLashof.zero())


def sQ2Q(m) -> DyerLashof:
    """ return sQ^m in terms of Q^I """
    if len(m) == 0:
        return DyerLashof.unit()
    elif len(m) == 1:
        return DyerLashof.gen(m[0])
    else:
        r = m[0]  # recursive
        d = sum(m[1:])
        return sum((DyerLashof.gen(r + i) * PQ(i, sQ2Q(m[1:]))
                    for i in range(max((d - r + 1) // 2, 0), d // 2 + 1)), DyerLashof.zero())


# ------------------ convert Q^I to sQ^J --------------------------------
def _my_indices(d: int, m: tuple):
    """ return an iterator of tuples for the Cartan formula for chiPsQ_mon """
    if len(m) == 0:
        return
    elif len(m) == 1:
        if m[0] >= d and binom_mod2(m[0], d):
            yield m[0] - d,
        return
    for i in range(min(d, m[0]) + 1):
        if binom_mod2(m[0], i):
            for t in _my_indices(d - i, m[1:]):
                yield (m[0] - i,) + t


def chiPsQ_mon(r: int, m: tuple) -> MyDyerLashof:
    """ return chiP^r sQ^m truncated by 0 """
    return MyDyerLashof(set(_my_indices(r, m)))


def chiPsQ(r: int, x: MyDyerLashof) -> MyDyerLashof:
    return sum((chiPsQ_mon(r, m) for m in x.data), MyDyerLashof.zero())


def Q2sQ(m) -> MyDyerLashof:
    """ return Q^m in terms of sQ^I """
    if len(m) == 0:
        return MyDyerLashof.unit()
    elif len(m) == 1:
        return MyDyerLashof.gen(m[0])
    else:
        r = m[0]  # recursive, complexity = length
        d = sum(m[1:])
        return sum((MyDyerLashof.gen(r + i) * chiPsQ(i, Q2sQ(m[1:]))
                    for i in range(0, d + 1)), MyDyerLashof.zero()).simplify()


# ----------------- other methods --------------------------------------
def simplify_sQ(I):
    """ make sQ^I a sum of admissible sQ^J """
    result = set()
    residue = sQ2Q(I)
    while residue:
        m = max(residue.data)
        m_sq = tuple((m[i] - sum(m[i+1:])) * (1 << i) for i in range(len(m)))
        residue -= sQ2Q(m_sq)
        result.add(m_sq)
    return result


# ----------------- test -------------------------------------
def ordered_indexing_set(I, k, max_jn=None):
    """ iterator of J such that deg(J)=k and J<=I"""
    if max_jn is None:
        max_jn = k
    if len(I) == 1:
        if k <= I[0] and k <= max_jn:
            yield k,
        return
    for jn in range((k + 2) // 3, min(I[-1], k, max_jn) + 1):
        for J in ordered_indexing_set(I[:-1], k - jn, jn):
            yield J + (jn,)


if __name__ == "__main__":
    pass
