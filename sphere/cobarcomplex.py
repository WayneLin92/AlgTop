"""This module implements the cobarcomplex of the dual Steenrod algebra."""
import itertools
from algebras.mymath import Vector, orderedpartition
import algebras.BaseAlgebras as BA
from algebras.operations import DualSteenrod
from algebras.mymath import two_expansion
# TODO: autocomplete the representing cycle


class CobarSteenrod(BA.AlgebraMod2):
    # ---------- AlgebraMod2 --------------
    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple):
        return mon1 + mon2

    @staticmethod
    def deg_mon(mon: tuple):
        return Vector((len(mon), sum(DualSteenrod.deg_mon(m) for m in mon)))

    @staticmethod
    def str_mon(mon: tuple):
        if len(mon) == 0:
            return "1"
        result = "["
        for i in reversed(range(len(mon))):
            result += CobarSteenrod.str_mon_xi(mon[i])  #
            if i > 0:
                result += "|"
        result += "]"
        return result

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

    # Methods -----------------
    @classmethod
    def tensor(cls, *args: DualSteenrod):
        """Return $a_1\\otimes\\cdots\\otimes a_n$."""
        multi_data = reversed(tuple(r.data for r in args))
        iter_prod = itertools.product(*multi_data)
        data = set(m for m in iter_prod)
        return cls(data)

    @staticmethod
    def str_mon_xi(m):
        result = ""
        for g, e in m:
            for s in two_expansion(e):
                result += f"\\xi_{{{s}{s + g}}}"
        if result == "":
            result = "1"
        return result

    @staticmethod
    def _xi(*args):
        """Return the monomial $\\prod\\xi_{t_i-s_i}^{2^{s_i}}$."""
        result = DualSteenrod.unit()
        for i in range(len(args) // 2):
            s, t = args[2*i], args[2*i+1]
            result *= DualSteenrod.gen(t - s, 2 ** s)
        for m in result.data:
            return m

    @classmethod
    def R(cls, i, j):
        R"""Return [\xi^{i}_{j-i}]"""
        return DualSteenrod.gen(j - i, 2 ** i)

    @classmethod
    def h(cls, i: int, S: tuple = ()):
        """Return the representing cycle for h_i(S)."""
        n = len(S) + 1
        seq = set(range(i, i + 2 * n))
        S = {i + s for s in S} | {i}
        T = seq - S
        assert len(S) + len(T) == 2 * n
        S, T = sorted(S), sorted(T)
        data = set()
        for T1 in itertools.permutations(T):
            if all(t - s > 0 for s, t in zip(S, T1)):
                m = tuple(cls._xi(s, t) for s, t in reversed(tuple(zip(S, T1))))
                data.add(m)

                for i, j in itertools.combinations(range(n), 2):
                    if S[j] < T1[i] < T1[j]:
                        m = tuple(cls._xi(S[i], S[j]) if ell == i else
                                  cls._xi(S[j], T1[i], S[j], T1[j]) if ell == j else
                                  cls._xi(S[ell], T1[ell]) for ell in reversed(range(n)))
                        data.add(m)
                    if S[j] < T1[i]:
                        for k in range(i + 1, j):
                            m = tuple(cls._xi(S[i], S[j]) if ell == i else
                                      cls._xi(S[j], T1[i], S[k], T1[k]) if ell == k else
                                      cls._xi(S[ell], T1[ell]) for ell in reversed(range(n)))
                            data.add(m)
                    if S[i] < T1[j] < T1[i]:
                        for k in range(i + 1, j + 1):
                            m = tuple(cls._xi(S[i], T1[j]) if ell == i else
                                      cls._xi(T1[j], T1[i], S[k], T1[k]) if ell == k else
                                      cls._xi(S[ell], T1[ell]) for ell in reversed(range(n)))
                            data.add(m)
        return cls(data)

    @classmethod
    def b(cls, i, j):
        """Return the representing cycle for b_{ij}."""
        data = {(cls._xi(i, j), cls._xi(i, j))}
        for k in range(i + 1, j):
            mon1 = (cls._xi(i, j, k, j), cls._xi(i, k))
            mon2 = (cls._xi(k, j), cls._xi(i, j, i, k))
            data ^= {mon1, mon2}
        return cls(data)

    @staticmethod
    def basis_mon(s, t, u):
        """Return a basis of degree s, t and weight <= u."""
        for d in orderedpartition(s, t):
            iters = (DualSteenrod.basis_mons(d[i]) for i in range(s))
            for mon in itertools.product(*iters):
                if CobarSteenrod.weight_mon(mon) <= u:
                    yield mon

    @staticmethod
    def weight_mon(mon: tuple):
        return sum(DualSteenrod.weight_mon(m) for m in mon)

    def weight(self):
        return max(map(self.weight_mon, self.data))

    def terms_topweight(self):
        tw = self.weight()
        data = {m for m in self.data if self.weight_mon(m) == tw}
        return type(self)(data)

    @staticmethod
    def is_simple(mon):
        for m in mon:
            if len(m) != 1 or bin(m[0][1]).count("1") != 1:
                return False
        return True

    def terms_simple(self):
        return type(self)({mon for mon in self.data if self.is_simple(mon)})

    def diff(self):
        """Return the differential."""
        data = set()
        for mon in self.data:
            for i in range(len(mon)):
                m = mon[i]
                coprod = DualSteenrod(m).coprod()
                for m1, m2 in coprod.data:
                    if m1 and m2:
                        data ^= {mon[:i] + (m1, m2) + mon[i+1:]}
        return type(self)(data)



