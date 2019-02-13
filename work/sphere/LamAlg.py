from algebras import BaseAlgebras as BC
from algebras.mymath import binom_mod2
import algebras.linalg as linalg


class Lambda_Algebra(BC.OperationsMod2):
    """ This is for the Lambda algebra of lambda^I at prime 2 """
    # -- OperationsMod2 --------------
    @staticmethod
    def str_mon(mon: tuple):
        str_result = ""
        for i in mon:
            if i >= 10:
                str_result += "\\lambda_{{{0}}}".format(i)
            else:
                str_result += "\\lambda_{0}".format(i)
        if mon == ():
            str_result = "1"
        return str_result

    @staticmethod
    def is_null(mon, degree=None):
        return False

    @staticmethod
    def is_admissible(mon):
        for i in range(len(mon) - 1):
            if 2 * mon[i] < mon[i + 1]:
                return False, i
        return True, None

    @staticmethod
    def adem(mon, i):
        return set((mon[:i] + (mon[i+1] - mon[i] - 1 - k, 2 * mon[i] + 1 + k) + mon[i+2:])
                   for k in range(mon[i+1] // 2 - mon[i])
                   if binom_mod2(mon[i+1] - 2 * (mon[i] + 1 + k), k) == 1)

    @classmethod
    def gen(cls, *n: int) -> "Lambda_Algebra":
        return cls(n).simplify()

    # methods ----------------------
    @staticmethod
    def basis_mons(length, deg, i_min=0):
        this_cls = Lambda_Algebra
        if length == 0:
            if deg == 0:
                yield ()
            return
        if length == 1:
            if deg - 1 >= i_min:
                yield (deg - 1,)
            return
        for i in range(i_min, deg):
            for t in this_cls.basis_mons(length - 1, deg - i - 1, (i + 1) // 2):
                yield t + (i,)

    @classmethod
    def basis(cls, length, deg):
        return (cls(m) for m in Lambda_Algebra.basis_mons(length, deg))

    def diff(self):
        data = set()
        for m in self.data:
            for i in range(len(m)):
                data ^= {m[:i] + (m[i] - j, j - 1) + m[i+1:]
                         for j in range(1, m[i] // 2 + 1)
                         if binom_mod2(m[i] - 2 * j, j)}
        return Lambda_Algebra(data).simplify()

    @classmethod
    def homology(cls, s, t):
        assert s > 0
        my_map1 = linalg.LinearMapKernelMod2()
        my_map2 = linalg.LinearMapKernelMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t))
        print("kernel dim:", my_map1.kernel.dim())
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s - 1, t))
        print("image: dim", my_map2.image.dim())
        print("quotient:")
        for r in my_map1.kernel.quotient(my_map2.image).basis(Lambda_Algebra):
            print(r)

