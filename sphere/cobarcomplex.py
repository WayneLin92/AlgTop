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
            result += CobarSteenrod.str_mon_R(mon[i])
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
    def str_mon_R(m):
        result = ""
        for g, e in m:
            for s in two_expansion(e):
                result += f"R_{{{s}{s + g}}}"
        if result == "":
            result = "1"
        return result

    @staticmethod
    def R_mon(*args):
        R"""Return the monomial $\prod R_{ij}$.

        args is i1, j1, i2, j2, ..."""
        result = DualSteenrod.unit()
        for i in range(len(args) // 2):
            s, t = args[2*i], args[2 * i + 1]
            result *= DualSteenrod.gen(t - s, 2 ** s)
        for m in result.data:
            return m

    @staticmethod
    def R(i, j):
        """Return $R_{ij}$"""
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
                m = tuple(cls.R_mon(s, t) for s, t in reversed(tuple(zip(S, T1))))
                data.add(m)
        return cls(data)

    @classmethod
    def b(cls, i, j):
        """Return the representing cycle for b_{ij}."""
        data = {(cls.R_mon(i, j), cls.R_mon(i, j))}
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

    def summands_topweight(self):
        tw = self.weight()
        data = {m for m in self.data if self.weight_mon(m) == tw}
        return type(self)(data)

    @staticmethod
    def is_simple(mon):
        for m in mon:
            if len(m) != 1 or bin(m[0][1]).count("1") != 1:
                return False
        return True

    def summands_simple(self):
        return type(self)({mon for mon in self.data if self.is_simple(mon)})

    @staticmethod
    def key_mon(mon):
        return [(DualSteenrod.is_gen_E0(m), m) for m in mon]

    @staticmethod
    def _get_i(mon):
        for i in range(len(mon)):
            if DualSteenrod.is_gen_E0(mon[i]):
                if i > 0 and mon[i] < mon[i - 1]:
                    return i - 1
            else:
                if i > 0 and all((g, 1 << s) < mon[i - 1][0] for g, e in mon[i] for s in two_expansion(e)):
                    return i - 1
                else:
                    break
        return None

    @classmethod
    def d0_inv_data(cls, data: set):
        """Find a cycle c such that $d_0c = data$."""
        data = data.copy()
        result = set()
        while data:
            mon = max(data, key=cls.key_mon)
            i = cls._get_i(mon)
            if i is None:
                print(cls(tuple(reversed(mon))))
                raise BA.MyValueError("Not d0 invertible")
                print("Not d0 invertible")
                break
            m_d0_inv = mon[:i] + (DualSteenrod.mul_mons(mon[i + 1], mon[i]),) + mon[i + 2:]
            assert m_d0_inv not in result
            result ^= {m_d0_inv}
            data ^= cls(m_d0_inv).d0().data
            BA.Monitor.print(len(data))
        return result

    def d0_inv(self):
        return type(self)(self.d0_inv_data(self.data))

    def chain_d3(self):
        """Add chains to make it a d3 chain."""
        b1 = self.diff()
        a1 = b1.d0_inv()
        b2 = (self + a1).diff().summands_topweight()
        a2 = b2.d0_inv()
        return self + a1 + a2

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

    def d0(self):
        """Return the d_0 differential."""
        data = set()
        for mon in self.data:
            for i in range(len(mon)):
                m = mon[i]
                coprod = DualSteenrod(m).coprod_E0()
                for m1, m2 in coprod.data:
                    if m1 and m2:
                        data ^= {mon[:i] + (m1, m2) + mon[i+1:]}
        return type(self)(data)
