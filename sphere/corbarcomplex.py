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
    def weight_mon(mon: tuple):
        pass

    @classmethod
    def tensor(cls, *args: DualSteenrod):
        multi_data = tuple(r.data for r in args)
        iter_prod = itertools.product(*multi_data)
        data = set(m for m in iter_prod)
        return cls(data)

    def diff(self):
        data = set()
        for mon in self.data:
            for i in range(len(mon)):
                m = mon[i]
                coprod = DualSteenrod(m).coprod()
                for m1, m2 in coprod.data:
                    if m1 and m2:
                        data ^= {mon[:i] + (m1, m2) + mon[i+1:]}
        return type(self)(data)
