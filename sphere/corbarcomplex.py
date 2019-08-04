import itertools
from algebras.mymath import Deg
import algebras.BaseAlgebras as BA
from algebras.operations import DualSteenrod


class CobarSteenrod(BA.AlgebraMod2):
    # ---------- AlgebraMod2 --------------
    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple):
        return mon1 + mon2

    @staticmethod
    def deg_mon(mon: tuple):
        return Deg((len(mon), sum(DualSteenrod.deg_mon(m) for m in mon)))

    @staticmethod
    def str_mon(mon: tuple):
        if len(mon) == 0:
            return "1"
        result = "["
        for i in range(len(mon)):
            result += DualSteenrod.str_mon(mon[i])
            if i < len(mon) - 1:
                result += "|"
        result += "]"
        return result

    # Methods -----------------
    @staticmethod
    def basis_mon(s, t, u):
        """Return a basis of degree s, t and weight u."""


    @staticmethod
    def weight_mon(mon: tuple):
        return sum(DualSteenrod.weight_mon(m) for m in mon)

    @classmethod
    def tensor(cls, *args: DualSteenrod):
        """Return $a_1\\otimes\\cdots\\otimes a_n$."""
        multi_data = tuple(r.data for r in args)
        iter_prod = itertools.product(*multi_data)
        data = set(m for m in iter_prod)
        return cls(data)

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

    @staticmethod
    def cup(x, y, j):
        """Return $x\\cup_j y$."""
        data = set()
        for alpha in x.data:
            for beta in y.data:
                p, q = len(alpha), len(beta)
                if j % 2 == 0:
                    i = j // 2
                    for r1 in itertools.combinations(range(p-i-1), i):
                        for s1 in itertools.combinations(range(1, q-i), i):
                            r = [r1[k] + k + 1 for k in range(i)]
                            s = [s1[k] + k + 1 for k in range(i)]
                            m = [alpha[k] for k in range(r[0])]
                            for ell in range(i):
                                pass



