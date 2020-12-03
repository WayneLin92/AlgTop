"""Steenrod operations and Dyer-Lashof operations"""
import math
from typing import Union, List
from algebras import BaseAlgebras as BA
from algebras.mymath import choose_mod2, tex_pow, tex_sub, add_tuple, sub_tuple
# todo: add operations at odd primes
# todo: use mul_data instead of mul


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
        result = "".join(tex_pow('Q', i) for i in mon)
        if result == "":
            result = "1"
        return result

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

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

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
    def gen(cls, key: tuple = ()):
        return cls(((key, 1),))

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

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
    @staticmethod
    def repr_mon(mon, clsname="Steenrod"):
        return f"{clsname}({mon})"

    def __mul__(self, other):
        if type(other) is DyerLashof:
            return AR.tensor(DyerLashof(()), self) * other
        else:
            return super().__mul__(other)

    @staticmethod
    def str_mon(mon: tuple):
        result = "".join(f"\\mathit{{Sq}}^{i}" if i < 10 else f"Sq^{{{i}}}" for i in mon)
        if result == "":
            result = "1"
        return result

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
                # noinspection PyTypeChecker
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

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

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

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

    # -- HopfAlgebra --------------
    @classmethod
    def coprod_gen(cls, n):
        data = {(cls._mon_gen(n - i, 1 << i), cls._mon_gen(i)) for i in range(n + 1)}
        return DualSteenrodT2(data)

    @classmethod
    def coprod_gen_E0(cls, n):
        R"""Coproduct of $\xi_n$ in $E^0A_*$"""
        data = {((), cls._mon_gen(n)), (cls._mon_gen(n), ())}
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
    def is_gen_E0(mon):
        return len(mon) == 1 and bin(mon[0][1]).count('1') == 1

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

    def coprod_E0(self):
        """Coproduct in $E_0A_*$"""
        result = self.type_T2().zero()
        for m in self.data:
            product = result.unit()
            for gen, exp in m:
                product = product * (self.coprod_gen_E0(gen) ** exp)
            result += product
        return result

    @staticmethod
    def type_dual():
        return Steenrod

    @staticmethod
    def type_T2():
        return DualSteenrodT2

    # methods ---------------------
    @staticmethod
    def weight_mon(mon: tuple):
        """Return the weight of the the monomial."""
        return sum((2 * g - 1) * bin(e).count('1') for g, e in mon)

    def actQ(self, s):
        pass

    def actSq(self, s):
        pass

    @staticmethod
    def _mon_gen(n, e=1):
        if n < 0:
            raise ValueError(f"n(={n}) should be nonnegative")
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
        return (cls(m) for m in cls.basis_mons(deg))


class DualSteenrodT2(BA.AlgebraT2Mod2):
    """ Tensor product of two DualSteenrod """
    type_c0 = DualSteenrod
    type_c1 = DualSteenrod

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass


class AR(BA.AlgebraT2Mod2):
    """
    The AR algebra
    self.data is a set of (Q^I, Sq^I)
    """
    type_c0 = DyerLashof
    type_c1 = Steenrod

    # -- AlgebraMod2 ----------

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

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


class SteenrodMilnor(BA.AlgebraMod2):
    """Steenrod algebra using the Milnor basis."""

    # -- AlgebraMod2 -----------------
    @staticmethod
    def mul_mons(mon1, mon2):
        result = set()
        for X in SteenrodMilnor._get_Xs(mon1, mon2):
            T = X[0][1:]
            for i, row in enumerate(X[1:]):
                T = add_tuple(T, (0,) * i + row)
            result ^= {T}
        return result

    @staticmethod
    def str_mon(mon) -> str:
        return f"\\mathit{{Sq}}({', '.join(map(str, mon))})" if mon else "1"

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

    @staticmethod
    def deg_mon(mon: tuple):
        return sum((e << i + 1) - e for i, e in enumerate(mon))

    # Methods -----------------------
    @classmethod
    def gen(cls, *args: int):
        return cls(tuple(args))

    @staticmethod
    def _get_rows(r, S, j_max=None, allow_trailing_zeros=False, B=()):
        R"""Return (x_1,...) such that $\sum 2^ix_i=r$ and $x_i<=s_{i-1}$ and not x_i & b_i"""
        if j_max is None:
            j_max = len(S)
        if j_max == 0:
            yield r,
            return
        s = S[j_max - 1] if len(S) >= j_max else 0
        x = min(s, r >> j_max)
        for x_j_max in reversed(range(x + 1)):
            if len(B) <= j_max or not x_j_max & B[j_max]:
                for X in SteenrodMilnor._get_rows(r - (x_j_max << j_max), S, j_max - 1,
                                                  allow_trailing_zeros or x_j_max, B):
                    yield X + (x_j_max,) if allow_trailing_zeros or x_j_max else X

    @staticmethod
    def _get_Xs(R, S, B=()):
        if not R:
            if all(not s & b for s, b in zip(S, B[1:])):
                yield (None,) + S,
            return
        for row in SteenrodMilnor._get_rows(R[-1], S, None, False, B):
            for X in SteenrodMilnor._get_Xs(R[:-1], sub_tuple(S, row[1:]), (0,) + add_tuple(B, row)):
                yield X + (row,)


class DualSteenrodDense(BA.AlgebraMod2):
    """Dual Steenrod algebra with dense data structure."""

    # -- AlgebraMod2 -------------
    @staticmethod
    def mul_mons(mon1, mon2):
        return add_tuple(mon1, mon2)

    @staticmethod
    def str_mon(mon) -> str:
        return "".join(tex_pow(tex_sub("\\xi", i + 1), e) for i, e in enumerate(mon) if e) if mon else "1"

    @staticmethod
    def deg_mon(mon: tuple):
        return sum(e * ((1 << i + 1) - 1) for i, e in enumerate(mon))

    @classmethod
    def gen(cls, n: int, e: int = 1):
        return cls(cls._mon_gen(n, e))

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass

    # -- HopfAlgebra --------------
    @classmethod
    def coprod_gen(cls, n):
        data = {(cls._mon_gen(n - i, 1 << i), cls._mon_gen(i)) for i in range(n + 1)}
        return cls.type_T2()(data)

    @classmethod
    def coprod_gen_E0(cls, n):
        R"""Coproduct of $\xi_n$ in $E^0A_*$"""
        data = {((), cls._mon_gen(n)), (cls._mon_gen(n), ())}
        return cls.type_T2()(data)

    @staticmethod
    def is_gen(mon): #???
        if sum(map(bool, mon)) == 1:
            return True, len(mon)
        else:
            return False, None

    @classmethod
    def is_gen_E0(cls, mon):
        return cls.is_gen(mon)[0] and bin(mon[-1]).count('1') == 1

    def coprod(self):
        type_T2 = self.type_T2()
        result = type_T2.zero_data()
        for m in self.data:
            product = type_T2.unit_data()
            for i, e in enumerate(m):
                product = type_T2.mul_data(product, (self.coprod_gen(i + 1) ** e).data)
            result ^= product
        return type_T2(result)

    def coprod_E0(self):
        type_T2 = self.type_T2()
        result = type_T2.zero_data()
        for m in self.data:
            product = type_T2.unit_data()
            for i, e in enumerate(m):
                product = type_T2.mul_data(product, (self.coprod_gen_E0(i + 1) ** e).data)
            result ^= product
        return type_T2(result)

    @staticmethod
    def type_T2():
        return DualSteenrodDenseT2

    # methods ---------------------
    @staticmethod
    def weight_mon(mon: tuple):
        """Return the weight of the the monomial."""
        return sum((2 * j + 1) * bin(e).count('1') for j, e in enumerate(mon))

    @staticmethod
    def _mon_gen(n, e=1):
        if n < 0:
            raise ValueError(f"n(={n}) should be nonnegative")
        return (0,) * (n - 1) + (e,) if n > 0 else ()

    @classmethod
    def basis_mons(cls, deg, index_max=None):
        if deg == 0:
            yield ()
            return
        if index_max is None:
            index_max = int(math.log2(deg + 1))
        if index_max == 1:
            yield (deg,) if deg > 0 else ()
        else:
            deg_index_max = (1 << index_max) - 1
            for i in range(deg // deg_index_max, -1, -1):
                for mon in cls.basis_mons(deg - i * deg_index_max, index_max - 1):
                    if i > 0:
                        yield mon + (0,) * (index_max - 1 - len(mon)) + (i,)
                    else:
                        yield mon

    @classmethod
    def basis(cls, deg):
        return (cls(m) for m in cls.basis_mons(deg))


class DualSteenrodDenseT2(BA.AlgebraT2Mod2):
    """ Tensor product of two DualSteenrod """
    type_c0 = DualSteenrodDense
    type_c1 = DualSteenrodDense

    @staticmethod
    def repr_mon(mon, clsname) -> str:
        pass


if __name__ == "__main__":
    xi = DualSteenrodDense.gen
    print(xi(3), xi(3).deg())
    for x in DualSteenrodDense.basis(10):
        print(x)
