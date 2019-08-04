from . import BaseAlgebras as BA
from .mymath import binom_mod2, cartanindices


# Classes -----------------------------------------------------
class HBOZ(BA.BasePolyMod2):
    """
    Homology of BO\\times Z.
    Should not be inherited.
    """
    # -- BasePolyMod2 ------------
    @classmethod
    def gen(cls, n: int):
        return cls(((n, 1),))

    @staticmethod
    def deg_gen(n: int) -> int:
        return n

    @staticmethod
    def str_gen(n: int) -> str:  # ######## b_0^0
        if n >= 10:
            return "b_{{{}}}".format(n)
        else:
            return "b_{}".format(n)

    # methods --------------------
    @staticmethod
    def sQ_gen(n, s):
        this_class = HBOZ
        result = set()
        for i in range(n + 1):
            if binom_mod2(s, i):
                result ^= {((s + i, 1), (n - i, 1))} if s + i < n - i \
                    else {((n - i, 1), (s + i, 1))} if s + i > n - i \
                    else {((s + i, 2),)}
        return this_class(result)

    @classmethod
    def sQ_gen_pow(cls, n, e, s) -> "HBOZ":
        if e == 1:
            return cls.sQ_gen(n, s)
        if e % 2:
            result = cls.zero()
            for i in range(0, s + 1, 2):  # recursive but fast
                result += cls.sQ_gen(n, s - i) * cls.sQ_gen_pow(n, (e - 1) // 2, i // 2).square()
            return result
        else:
            return cls.sQ_gen_pow(n, e // 2, s // 2).square() if s % 2 == 0 else cls.zero()

    def sQ(self, s):
        result = self.zero()
        for m in self.data:
            n = len(m)
            for t in cartanindices(n, s):
                product = self.unit()
                for i in range(n):
                    product *= self.sQ_gen_pow(m[i][0], m[i][1], t[i])
                result += product
        return result
