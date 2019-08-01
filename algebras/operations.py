""" Provides the classes for the Dyer-Lashof operations
    and the Steenrod operations and its dual
"""
from . import BaseAlgebras as BA
from .mymath import choose_mod2
import math
from typing import Union, List


# Classes -----------------------------------------------
class DyerLashof(BA.OperationsMod2):
    """ This is for the Dyer_Lashof algebra of Q^I at prime 2 on the space level """
    # -- OperationsMod2 --------------
    def __mul__(self, other):
        if type(other) is Steenrod:
            return AR.tensor(self, other)
        else:
            return super().__mul__(other)

    @staticmethod
    def str_mon(mon: tuple):
        str_result = ""
        for i in mon:
            if i >= 10:
                str_result += "Q^{{{0}}}".format(i)
            else:
                str_result += "Q^{0}".format(i)
        if mon == ():
            str_result = "1"
        return str_result

    @staticmethod
    def is_null(mon, degree=None):
        partial_sum = degree if degree is not None else 0
        if mon == ():
            return False
        if mon[-1] < partial_sum:
            return True
        for i in range(len(mon) - 1, 0, -1):
            partial_sum += mon[i]
            if mon[i - 1] < partial_sum:
                return True
        return False

    @staticmethod
    def is_admissible(mon):
        for i in range(len(mon) - 1):
            if mon[i] > 2 * mon[i + 1]:
                return False, i
        return True, None

    @staticmethod
    def adem(mon, index):
        return set((mon[:index] + (mon[index] + mon[index + 1] - j, j) + mon[index + 2:])
                   for j in range((mon[index] + 1) // 2, (mon[index] + mon[index + 1]) // 2 + 1)
                   if choose_mod2(j - mon[index + 1] - 1, 2 * j - mon[index]) == 1)

    # methods -------------
    @classmethod
    def gen(cls, *n: int):
        return cls(n).simplify()

    @staticmethod
    def basis_mons(deg, i_next=None, i_sum=0):
        this_cls = DyerLashof
        if deg == 0:
            yield ()
        elif i_next is None:
            for j in range(i_sum + 1, deg + 1):
                for t in this_cls.basis_mons(deg - j, j, j + i_sum):
                    yield t + (j,)
        else:
            for j in range(i_sum + 1, min(2 * i_next + 1, deg + 1)):
                for t in this_cls.basis_mons(deg - j, j, j + i_sum):
                    yield t + (j,)

    @classmethod
    def basis(cls, deg):
        return (cls(m) for m in cls.basis_mons(deg))

    @staticmethod
    def excess(m):
        """ Return the excess of the monomial """
        if m == ():
            return 0
        else:
            return m[0] - sum(m[1:])


class DyerLashofX(BA.BasePolyMod2):
    """
    This is for Dyer_lashof operations acting on a fix element
    Data is a set of tuples
    Each tuple is actually a dict
    """
    x_deg = 0

    # -- BasePolyMod2 ---------------
    @classmethod
    def deg_gen(cls, key: tuple):
        return DyerLashof.deg_mon(key) + cls.x_deg

    @staticmethod
    def str_gen(key: tuple):
        if key != ():
            return "{}x".format(DyerLashof(key))
        else:
            return "x"

    @classmethod
    def str_mon(cls, mon):
        result = ""
        for gen, exp in mon:
            if exp >= 10 or exp < 0:
                result += "({})^{{{}}}".format(cls.str_gen(gen), exp)
            elif exp > 1:
                result += "({})^{}".format(cls.str_gen(gen), exp)
            elif exp == 1:
                result += cls.str_gen(gen)
        if result == "":
            result = "1"
        return result

    @classmethod
    def gen(cls, key: tuple = ()):
        return cls(((key, 1),))

    # methods ----------------
    def __rmul__(self, other):
        if type(other) is type(self):
            return super().__mul__(other)
        elif type(other) is DyerLashof:
            return sum((self.actQ(m) for m in other.data), self.zero())
        else:
            return NotImplemented

    @classmethod
    def set_x_deg(cls, x_deg):
        cls.x_deg = x_deg

    @staticmethod
    def sqrt(m):
        """
        Return the square root of the monomial
        Return None if it is not a square~
        """
        for item in m:
            if item[1] % 2 == 1:
                return None
        return tuple((item[0], item[1] // 2) for item in m)

    @staticmethod
    def decompose(m):
        """
        Decompose a monomial into two monomials
        m1 is a generator or a power of a generator
        m must be decomposable
        """
        return (m[0:1], m[1:]) if len(m) > 1 else (((m[0][0], 1),), ((m[0][0], m[0][1] - 1),))

    def actQ(self, s):
        if type(s) is tuple:
            result = self
            for i in range(len(s) - 1, -1, -1):
                result = result.actQ(s[i])
            return result

        result = self.zero()
        for m in self.data:
            if len(m) == 1 and m[0][1] == 1:
                result += DyerLashofX((((s,) + m[0][0], 1),))
            elif len(m) == 0:
                if s == 0:
                    result += self.unit()
            else:
                mr = self.sqrt(m)
                if mr is not None:
                    if s % 2 == 0:
                        result += DyerLashofX(mr).actQ(s // 2).square()
                else:
                    m1, m2 = self.decompose(m)
                    result = sum((DyerLashofX(m1).actQ(i) * DyerLashofX(m2).actQ(s - i)
                                  for i in range(self.deg_mon(m1), s - self.deg_mon(m2) + 1)), result)
        return result.simplify()

    def simplify(self):
        result = self.zero()
        for mon in self.data:
            pro = self.unit()
            for k, v in mon:
                fk = DyerLashof(k).simplify(self.x_deg)
                pro = pro * DyerLashofX(set((self.simplify_square(m),) for m in fk.data) - {(None,)}) ** v
            result += pro
        return result

    def simplify_square(self, mon):
        """
        simplify monomial m using Q^q x = x^2 if deg(x) = q
        """
        partial_sum = self.x_deg
        indices = []
        for i in range(len(mon) - 1, -1, -1):
            if mon[i] == partial_sum:
                indices.append(i)
            partial_sum += mon[i]
        e = 1
        m1 = []
        for i in range(len(mon) - 1, -1, -1):
            if i in indices:
                e *= 2
            else:
                if mon[i] % e == 0:
                    m1.append(mon[i] // e)
                else:
                    return None
        return tuple(m1[::-1]), e


class Steenrod(BA.HopfAlgWithDualMod2, BA.OperationsMod2):
    """ This is for the Steenrod algebra of Sq^I at prime 2 on the space level """
    _chi_sq = []  # type: List["Steenrod"]

    # -- OperationsMod2 --------------
    def __mul__(self, other):
        if type(other) is DyerLashof:
            return AR.tensor(DyerLashof(()), self) * other
        else:
            return super().__mul__(other)

    @staticmethod
    def str_mon(mon: tuple):
        str_result = ""
        for i in mon:
            if i >= 10:
                str_result += "Sq^{{{0}}}".format(i)
            else:
                str_result += "Sq^{0}".format(i)
        if mon == ():
            str_result = "1"
        return str_result

    @staticmethod
    def is_null(mon, degree=None):
        """ determine if mon is zero for degree reason """
        if degree is None:  # todo: create a cache for which Sq^iSq^j=0
            return False
        partial_sum = degree
        if mon == ():
            return False
        if mon[len(mon) - 1] > partial_sum:
            return True
        for i in range(len(mon) - 1, 0, -1):
            partial_sum += mon[i]
            if mon[i - 1] > partial_sum:
                return True
        return False

    @staticmethod
    def is_admissible(mon):
        """
        Determine if mon is admissible and return also the
        place to apply the adem relations
        """
        for i in range(len(mon) - 1):
            if mon[i] < 2 * mon[i + 1]:
                return False, i
        return True, None

    @staticmethod
    def adem(mon, index):
        # BA.Monitor.count(mon)
        return set((mon[:index] + ((mon[index] + mon[index + 1] - k, k) if k > 0 else (mon[index] + mon[index + 1],)) +
                    mon[index + 2:]) for k in range(mon[index] // 2 + 1)
                   if choose_mod2(mon[index + 1] - k - 1, mon[index] - 2 * k) == 1)

    @classmethod
    def gen(cls, *n: int):
        m = tuple(i for i in n if i > 0)
        return cls(m).simplify()

    # -- HopfAlgebra --------------
    @classmethod
    def coprod_gen(cls, n):
        data = {(cls._mon_gen(i), cls._mon_gen(n - i)) for i in range(n + 1)}
        return SteenrodT2(data)

    @staticmethod
    def complexity(mon):
        return len(mon)

    @staticmethod
    def sqrt(mon):
        return None

    @staticmethod
    def virschiebung(mon):
        """
        Return the Virschiebung of the monomial of Steenrod.
        Return None if it is zero.
        """
        for i in mon:
            if i % 2 == 1:
                return None
        return tuple(i // 2 for i in mon)

    @staticmethod
    def is_gen(mon):
        if len(mon) == 1:
            return True, mon[0]
        else:
            return False, None

    @staticmethod
    def decompose(mon):
        return (mon[0],), mon[1:]

    @staticmethod
    def pair_gen(n, mon_dual):
        return 1 if len(mon_dual) == 1 and mon_dual[0][0] == 1 else 0

    def coprod(self):
        """ coproduct of Steenrod algebra~ """
        result = self.type_T2().zero()
        for m in self.data:
            product = result.unit()
            for i in m:
                product = product * self.coprod_gen(i)
            result += product
        return result

    @staticmethod
    def type_dual(): return DualSteenrod

    @staticmethod
    def type_T2(): return SteenrodT2

    # methods ----------------------
    @staticmethod
    def _mon_gen(n):
        return (n,) if n > 0 else ()

    @staticmethod
    def excess(m):
        """ Return the excess of the monomial """
        if m == ():
            return 0
        else:
            return m[0] - sum(m[1:])

    @staticmethod
    def conj_gen(n):
        cls = Steenrod
        """ Return chi Sq^n """
        if len(cls._chi_sq) > n:
            return cls._chi_sq[n]
        else:
            if len(cls._chi_sq) == 0:
                cls._chi_sq.append(cls.unit())
            for i in range(len(cls._chi_sq), n + 1):
                sq_i = sum((cls.gen(j) * cls._chi_sq[i - j] for j in range(1, i + 1)), cls.zero())
                cls._chi_sq.append(sq_i)
            return cls._chi_sq[n]

    @staticmethod
    def basis_mons(deg, i_min=1):
        this_cls = Steenrod
        if deg == 0:
            yield ()
        else:
            for i in range(i_min, deg + 1):
                for t in this_cls.basis_mons(deg - i, 2 * i):
                    yield t + (i,)

    @classmethod
    def basis(cls, deg):
        return (cls(m) for m in Steenrod.basis_mons(deg))


class SteenrodT2(BA.AlgebraT2Mod2):  # todo: T2 for odd primes
    """ Tensor product of two Steenrod algebras """
    type_c0 = Steenrod
    type_c1 = Steenrod

    def __mul__(self, other):
        return super().__mul__(other).simplify()

    def simplify(self, degree=None):  # #####################
        """ Simplify the tensor product of this non-commutative algebra """
        return sum((self.tensor(Steenrod(m1).simplify(degree), Steenrod(m2).simplify(degree))
                    for m1, m2 in self.data), self.zero())


class DualSteenrod(BA.HopfAlgWithDualMod2, BA.BasePolyMod2):
    """ This is for the dual Steenrod algebra of xi_n~ """
    # -- BasePolyMod2 -------------
    @classmethod
    def gen(cls, n: int, e: int = 1):
        return cls(cls._mon_gen(n, e))

    @staticmethod
    def deg_gen(n: int):
        return (1 << n) - 1

    @staticmethod
    def str_gen(n: int):
        if n >= 10:
            return "\\xi_{{{}}}".format(n)
        else:
            return "\\xi_{}".format(n)

    # -- HopfAlgebra --------------
    @classmethod
    def coprod_gen(cls, n):
        data = {(cls._mon_gen(n - i, 1 << i), cls._mon_gen(i)) for i in range(n + 1)}
        return DualSteenrodT2(data)

    @staticmethod
    def complexity(mon):
        return (sum(m[1] for m in mon) - 1) * 2 + 1

    @staticmethod
    def sqrt(mon):
        for m in mon:
            if m[1] % 2 == 1:
                return None
        return tuple((m[0], m[1] // 2) for m in mon)

    @staticmethod
    def virschiebung(mon: tuple):
        return None

    @staticmethod
    def is_gen(mon):
        if len(mon) == 1 and mon[0][1] == 1:
            return True, mon[0][0]
        else:
            return False, None

    @staticmethod
    def decompose(mon):
        if len(mon) == 1:
            return ((mon[0][0], 1),), ((mon[0][0], mon[0][1] - 1),)
        else:
            return (mon[0],), mon[1:]

    @staticmethod
    def pair_gen(n, mon_dual):
        return 1 if mon_dual == tuple(1 << i for i in range(n - 1, -1, -1)) else 0

    def coprod(self):
        result = self.type_T2().zero()
        for m in self.data:
            product = result.unit()
            for gen, exp in m:
                product = product * (self.coprod_gen(gen) ** exp)
            result += product
        return result

    @staticmethod
    def type_dual():
        return Steenrod

    @staticmethod
    def type_T2():
        return DualSteenrodT2

    # methods ---------------------
    def actQ(self, s):
        pass

    def actSq(self, s):
        pass

    @staticmethod
    def _mon_gen(n, e=1):
        if n < 0:
            raise IndexError("n(={}) should be nonnegative".format(n))
        return ((n, e),) if n > 0 else ()

    @staticmethod
    def basis_mons(deg, index_max=None):
        this_cls = DualSteenrod
        if deg == 0:
            yield ()
            return
        if index_max is None:
            index_max = int(math.log2(deg + 1))
        if index_max == 1:
            yield ((1, deg),) if deg > 0 else ()
        else:
            for i in range(deg // this_cls.deg_gen(index_max), -1, -1):
                for mon in this_cls.basis_mons(deg - i * this_cls.deg_gen(index_max), index_max - 1):
                    if i > 0:
                        yield mon + ((index_max, i),)
                    else:
                        yield mon

    @classmethod
    def basis(cls, deg):
        return (cls(m) for m in DualSteenrod.basis_mons(deg))


class DualSteenrodT2(BA.AlgebraT2Mod2):
    """ Tensor product of two DualSteenrod """
    type_c0 = DualSteenrod
    type_c1 = DualSteenrod


class AR(BA.AlgebraT2Mod2):
    """
    The AR algebra
    self.data is a set of (Q^I, Sq^I)
    """
    type_c0 = DyerLashof
    type_c1 = Steenrod

    # -- AlgebraMod2 ----------
    def __init__(self, data: Union[DyerLashof, Steenrod, set, tuple]):
        if type(data) is self.type_c0:
            self.data = set((m, ()) for m in data.data)
        if type(data) is self.type_c1:
            self.data = set(((), m) for m in data.data)
        else:
            super().__init__(data)

    def __mul__(self, other):
        if type(other) is Steenrod:
            return self * AR.tensor(DyerLashof(()), other)
        elif type(other) is DyerLashof:
            return self * AR.tensor(other, Steenrod(()))
        else:
            return BA.AlgebraMod2.__mul__(self, other)

    # -- AlgebraT2Mod2 -----------
    def mul_mons(self, mon1, mon2):
        if len(mon1[1]) == 0:
            return AR.tensor(self.type_c0(mon1[0]) * self.type_c0(mon2[0]), self.type_c1(mon2[1])).data
        elif len(mon2[0]) == 0:
            # Note that Sq^I is the dual action
            return AR.tensor(self.type_c0(mon1[0]), self.type_c1(mon2[1]) * self.type_c1(mon1[1])).data
        else:
            return (self.type_c0(mon1[0]) * AR(mon2).actSq(mon1[1])).data

    def str_mon(self, mon):
        l0 = self.type_c0.str_mon(mon[0])
        l1 = self.type_c1.str_mon(mon[1])
        if l0 == "1":
            return l1
        elif l1 == "1":
            return l0
        else:
            return l0 + l1

    # methods
    def __rmul__(self, other):
        if type(other) is Steenrod:
            return AR.tensor(DyerLashof.unit(), other) * self
        elif type(other) is DyerLashof:
            return AR.tensor(other, Steenrod.unit()) * self
        else:
            return NotImplemented

    def actSq(self, r):
        if type(r) is tuple:
            result = self
            for i in range(len(r)):  # Note that Sq^I is the dual action
                result = result.actSq(r[i])
            return result

        result = self.zero()
        for m in self.data:
            if m[0] == ():
                result += AR.tensor(DyerLashof.unit(), Steenrod(m[1]) * Steenrod.gen(r))
            else:
                s = m[0][0]
                if s >= r:  # recursive
                    result = sum((DyerLashof.gen(s - r + i) * AR((m[0][1:], m[1])).actSq(i)
                                  for i in range(max(0, r - (s + 1) // 2), r // 2 + 1)  # ####
                                  if choose_mod2(s - r, r - 2 * i) == 1), result)
        return result

    def truncate(self, deg_sq_max):
        """ truncate by the deg of Steenrod part """
        return AR(set(mon for mon in self.data if self.deg_mon(mon)[1] <= deg_sq_max))

    @staticmethod
    def sQ(r, k_max):
        return sum((DyerLashof.gen(r + k) * Steenrod.gen(k) for k in range(k_max + 1)), AR.zero())


# 712, 627, 621, 591, 583, 575, 599
# todo: add operations at odd primes
