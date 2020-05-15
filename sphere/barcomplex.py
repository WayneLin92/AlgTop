"""This module implements the barcomplex of the Steenrod algebra."""
from itertools import product, combinations
from algebras.mymath import Vector, orderedpartition
import algebras.BaseAlgebras as BA
from algebras.operations import SteenrodMilnor


class BarSteenrod(BA.AlgebraMod2):
    # ---------- AlgebraMod2 --------------
    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple):
        mon = mon1 + mon2
        result = set()
        n, n1, n2 = len(mon), len(mon1), len(mon2)
        for p in combinations(range(n), n1):
            pi = p + tuple(i for i in range(n) if i not in p)
            pi_inv = [pi.index(i) for i in range(n)]
            mon_pi = tuple(mon[pi_inv[i]] for i in range(n))
            result ^= {mon_pi}
        return result

    @staticmethod
    def deg_mon(mon: tuple):
        return Vector((len(mon), sum(SteenrodMilnor.deg_mon(m) for m in mon)))

    @staticmethod
    def str_mon(mon: tuple):
        if len(mon) == 0:
            return "1"
        result = R"\{"
        result += "|".join(map(SteenrodMilnor.str_mon, mon))
        result += R"\}"
        return result

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

    # Methods -----------------
    @classmethod
    def tensor(cls, *args: SteenrodMilnor):
        R"""Return $a_1\otimes\cdots\otimes a_n$."""
        multi_data = [r.data for r in args]
        iter_prod = product(*multi_data)
        data = set(m for m in iter_prod)
        return cls(data)

    @staticmethod
    def str_mon_sq(mon):
        pass

    @staticmethod
    def P(i, j):
        """Return P^i_{j-i}$."""
        return (0,) * (j - i - 1) + (i,)

    @classmethod
    def h(cls, i: int, S: tuple = ()):
        """Return the image of h_i(S) via embedding."""
        n = len(S) + 1
        seq = set(range(i, i + 2 * n))
        S = {i + s for s in S} | {i}
        T = seq - S
        assert len(S) + len(T) == 2 * n
        S, T = sorted(S), sorted(T)
        data = set()
        pass
        return cls(data)

    @staticmethod
    def basis_mon(s, t, u):
        """Return a basis of degree s, t and weight <= u."""
        for d in orderedpartition(s, t):
            iters = (SteenrodMilnor.basis_mons(d[i]) for i in range(s))
            for mon in product(*iters):
                if SteenrodMilnor.weight_mon(mon) <= u:
                    yield mon

    @staticmethod
    def weight_mon(mon: tuple):
        return sum(SteenrodMilnor.weight_mon(m) for m in mon)

    def weight(self):
        return max(map(self.weight_mon, self.data))

    def terms_topweight(self):
        tw = self.weight()
        data = {m for m in self.data if self.weight_mon(m) == tw}
        return type(self)(data)

    def diff(self):
        """Return the differential."""
        data = set()
        for mon in self.data:
            for i in range(len(mon) - 1):
                prod = SteenrodMilnor(mon[i]) * SteenrodMilnor(mon[i + 1])
                for m in prod.data:
                    data ^= {mon[:i] + (m,) + mon[i+2:]}
        return type(self)(data)
