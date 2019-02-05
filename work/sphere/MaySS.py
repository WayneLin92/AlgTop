import operator
import itertools
import pickle
from typing import Iterable
from algebras import BaseClasses as BC, linalg, mymath


class Signature(tuple):
    """A subclass of tuple modeling functions on nonnegative numbers"""
    def __new__(cls, iterable) -> "Signature":
        # noinspection PyTypeChecker
        return tuple.__new__(cls, iterable)

    def __add__(self, other):
        """Element-wise addition."""
        if len(self) < len(other):
            return Signature(itertools.chain(map(operator.add, self, other), other[len(self):]))
        else:
            return Signature(itertools.chain(map(operator.add, self, other), self[len(other):]))
        
    def __sub__(self, other):
        """Element-wise subtraction."""
        if len(self) < len(other):
            return Signature(itertools.chain(map(operator.sub, self, other), (-i for i in other[len(self):])))
        else:
            return Signature(itertools.chain(map(operator.sub, self, other), self[len(other):]))

    def __radd__(self, other):
        """This is implemented for supporting sum()."""
        return self.__add__(other) if other is not 0 else self

    def __mul__(self, other: int):
        """Broadcast multiplication."""
        return Signature(map(lambda x: x * other, self))

    def __bool__(self):
        """Return if any element in self is nonzero."""
        return any(self)

    def deg(self):
        """Return the degree of the monomial with self as the signature."""
        return sum(fn << n for n, fn in enumerate(self))

    def accumulate(self):
        """Return the sequence of the partial sum."""
        return Signature(itertools.accumulate(self[:-1]))

    def diff(self):
        """Return the sequence of self[i] - self[i-1]."""
        return Signature((self[0],) + tuple(self[i] - self[i-1] for i in range(1, len(self))) + (-self[-1],))


class MaySS(BC.BasePolyMod2):
    """ This is for the MSS starting from E_1 """
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *key: int):
        assert len(key) == 2 and key[0] > 0 <= key[1]
        return cls(((cls.ij2deg(key), 1),))

    @staticmethod
    def deg_gen(n: int) -> int:
        return n

    @staticmethod
    def str_gen(n: int) -> str:
        return "h_{{{}，{}}}".format(*MaySS.deg2ij(n))

    # methods -----------------
    @staticmethod
    def ij2deg(key: tuple) -> int:
        return (1 << key[0]) - 1 << key[1]

    @staticmethod
    def deg2ij(n: int) -> tuple:
        i = bin(n).count('1')
        j = n.bit_length() - i
        return i, j

    @staticmethod
    def basis_mons(length, deg, may_filtr, ij_max=None):
        """ return monomials h_ij^k...h_index_max^e with given length, deg, may_filtr """
        cls = MaySS
        if length == 0 or deg == 0 or may_filtr == 0:
            if length == deg == may_filtr:
                yield ()
            return
        if ij_max is None:
            bound = deg - length + 1
            if bound <= 0:
                return
            sum_ij = bound.bit_length()
            j = 0
            while (1 << sum_ij) - (1 << j) > bound:
                j += 1
            ij_max = (sum_ij - j, j)
        if ij_max == (1, 0):
            if length == deg == may_filtr:
                yield ((MaySS.ij2deg(ij_max), deg),)
            return
        for e in range(min(length, deg // MaySS.ij2deg(ij_max), may_filtr // ij_max[0]), -1, -1):
            # print(length, deg, may_filtr, "{}^{}".format(ij_max, e))
            index_next = (ij_max[0] - 1, ij_max[1] + 1) if ij_max[0] > 1 else (sum(ij_max) - 1, 0)
            for mon in cls.basis_mons(length - e, deg - e * MaySS.ij2deg(ij_max),
                                      may_filtr - e * ij_max[0], index_next):
                if e > 0:
                    yield mon + ((MaySS.ij2deg(ij_max), e),)
                else:
                    yield mon

    @classmethod
    def basis(cls, length, deg, may_filtr):
        return (cls(m) for m in MaySS.basis_mons(length, deg, may_filtr))

    def diff(self):
        """ return the coboundary of the cochain """
        result = self.zero()
        for m in self.data:
            for ind in range(len(m)):
                g, e = m[ind]
                if e % 2:
                    h = self.gen
                    i, j = MaySS.deg2ij(g)
                    dh_k = sum((h(i - l, j + l) * h(l, j) for l in range(1, i)), self.zero())
                    if e - 1:
                        result += MaySS(m[:ind] + m[ind + 1:]) * MaySS(((g, e - 1),)) * dh_k
                    else:
                        result += MaySS(m[:ind] + m[ind + 1:]) * dh_k
        return result

    @classmethod
    def homology(cls, s, t, u):
        my_map1 = linalg.LinearMapKernelMod2()
        my_map2 = linalg.LinearMapKernelMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        print("kernel dim:", my_map1.kernel.dim())
        # for r in my_map1.kernel.basis(MaySS):
        #     print(r)
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s - 1, t, u))
        print("image: dim", my_map2.image().dim())
        # for r in my_map2.image().basis(MaySS):
        #     print(r)
        print("quotient:")
        for r in my_map1.kernel.quotient(my_map2.image()).basis(MaySS):
            print(r)


class DualMaySS(BC.AlgebraMod2):
    """ This is the dual of the May spectral sequence """
    _maps = None
    _homology = None
    _loaded = False

    @classmethod
    def load(cls, s_max, t_max, u_max):
        if not cls._loaded:
            try:
                with open("MaySS_DualMaySS.pickle", "rb") as f:
                    cls._maps, cls._homology = pickle.load(f)
                    cls._loaded = True
            except FileNotFoundError:
                cls._maps = {}
                cls._homology = {}

        for s in range(s_max + 1, -1, -1):
            for t in range(s, t_max + 1):
                for u in range(s, u_max + 1):
                    if (s, t, u) not in cls._maps:
                        cls._maps[(s, t, u)] = linalg.LinearMapKernelMod2()
                        cls._maps[(s, t, u)].add_maps((r, r.diff()) for r in cls.basis(s, t, u))
                        if s <= s_max:
                            cycles = cls._maps[(s, t, u)].kernel
                            if (s + 1, t, u) in cls._maps:
                                boundaries = cls._maps[(s + 1, t, u)].image()
                                cls._homology[(s, t, u)] = list(cycles.quotient(boundaries).simplify().basis(set))
                            else:
                                cls._homology[(s, t, u)] = list(cycles.simplify().basis(set))
                            if not cls._homology[(s, t, u)]:
                                del cls._homology[(s, t, u)]

    @classmethod
    def save(cls):
        if cls._loaded:
            with open("MaySS_DualMaySS.pickle", "wb") as f:
                pickle.dump((cls._maps, cls._homology), f)

    # ----- AlgebraMod2 ----------------
    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple):
        result = dict(mon1)
        for ij, r in mon2:
            if ij in result:
                if result[ij] & r:
                    return set()
                else:
                    result[ij] += r
            else:
                result[ij] = r
        return tuple(sorted(result.items()))

    @classmethod
    def gen(cls, i, j, r=1) -> "DualMaySS":
        """Return gamma_r P^i_j."""
        return cls((((i, j), r),)) if r else cls(())

    @staticmethod
    def str_gen(item: tuple) -> str:
        ij, r = item
        return "\\gamma_{}(\\bar{{P}}^{}_{})".format(*map(mymath.tex_index, (r, *ij)))

    @classmethod
    def str_mon(cls, mon: tuple):
        mon = cls.lexicographic(mon)
        result = "".join(map(cls.str_gen, mon))
        return result if result else "1"

    @staticmethod
    def deg_gen(item: tuple):
        ij, r = item
        return mymath.Deg((r, ((1 << ij[1]) - 1 << ij[0]) * r, ij[1] * r))

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(cls.deg_gen, mon), mymath.Deg((0, 0, 0)))

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (self.deg_mon(m), m))

    # ----- Signature ------------------
    @staticmethod
    def sig_mon(mon: tuple):
        """Return the signature function of a monomial."""
        len_sig = max(k[0] + k[1] for k, r in mon)
        f = [0] * (len_sig + 1)
        for k, r in mon:
            f[k[0]] += r
            f[k[0] + k[1]] -= r
        return Signature(itertools.accumulate(f[:-1]))

    def sig(self):
        """Return the signature of self assuming that it is homogeneous in signature."""
        for m in self.data:
            return self.sig_mon(m)

    @staticmethod
    def _fixed_i_sig_mon(sig: Signature, i, j_max=None):
        """Return a monomial of gamma_r P^i_j for fixed i"""
        if sig[i] == 0:
            yield ()
            return
        if j_max is None:
            j_max = 1
            while i + j_max < len(sig) and sig[i + j_max] > 0:
                j_max += 1
        for j in range(j_max, 1, -1):
            r_max = min(sig[i:i + j])
            for r in range(1, r_max + 1):
                for m in DualMaySS._fixed_i_sig_mon(sig - ((0,) * i + (r,) * j), i, j - 1):
                    yield m + (((i, j), r),)
        j, r = 1, sig[i]
        if r > 0:
            yield (((i, j), r),)

    @staticmethod
    def basis_sig_mon(sig: Signature):
        if not sig:
            yield ()
            return
        i = 0
        while not sig[i]:
            i += 1
        for m in DualMaySS._fixed_i_sig_mon(sig, i):
            sig_m = DualMaySS.sig_mon(m)
            for m1 in DualMaySS.basis_sig_mon(sig - sig_m):
                yield m + m1

    @classmethod
    def basis_sig(cls, sig: Signature):
        return map(cls, cls.basis_sig_mon(sig))

    @classmethod
    def homology_signature(cls, sig: Signature):
        s_min = sum(i for i in sig.diff() if i > 0)
        s_max = sum(sig)
        lin_maps = {}
        for s in range(s_min, s_max + 2):
            lin_maps[s] = linalg.LinearMapKernelMod2()
        for x in cls.basis_sig(sig):
            s = x.deg()[0]
            lin_maps[s].add_maps(((x, x.diff()),))
        homology = {}
        for s in range(s_min, s_max + 1):
            print(f"{s}:")
            homology[s] = list(lin_maps[s].kernel.quotient(lin_maps[s+1].image()).simplify().basis(DualMaySS))
            for x in homology[s]:
                print(x)
        return homology

    # methods -----------------
    @staticmethod
    def ij2deg(ij: tuple) -> int:
        return (1 << ij[1]) - 1 << ij[0]

    @staticmethod
    def basis_mons(deg_s, deg_t, deg_u, ij_max=None):
        """ return monomials h_ij^k...h_index_max^e with given length, deg, may_filtr """
        if deg_s == 0 or deg_t == 0 or deg_u == 0:
            if deg_s == deg_t == deg_u:
                yield ()
            return
        if ij_max is None:
            bound = deg_t - deg_s + 1
            if bound <= 0:
                return
            sum_ij = bound.bit_length()
            i = 0
            while (1 << sum_ij) - (1 << i) > bound:
                i += 1
            ij_max = (i, sum_ij - i)
        if ij_max == (0, 1):
            if deg_s == deg_t == deg_u:
                yield (((0, 1), deg_t),)
            return
        for e in range(min(deg_s, deg_t // DualMaySS.ij2deg(ij_max), deg_u // ij_max[1]), -1, -1):
            index_next = (ij_max[0] + 1, ij_max[1] - 1) if ij_max[1] > 1 else (0, sum(ij_max) - 1)
            for mon in DualMaySS.basis_mons(deg_s - e, deg_t - e * DualMaySS.ij2deg(ij_max),
                                            deg_u - e * ij_max[1], index_next):
                if e > 0:
                    yield tuple(sorted(mon + ((ij_max, e),)))
                else:
                    yield mon

    @classmethod
    def basis(cls, deg_s, deg_t, deg_u):
        return map(cls, DualMaySS.basis_mons(deg_s, deg_t, deg_u))

    @staticmethod
    def deg_s_mon(mon):
        return sum(r for k, r in mon)

    @staticmethod
    def coprod_gen(item):
        deg, r = item
        yield (), (deg, r)
        for k in range(1, r):
            yield (deg, k), (deg, r-k)
        yield (deg, r), ()

    def coprod(self):
        data = set()
        for m in self.data:
            for comb in itertools.product(*map(self.coprod_gen, m)):
                items1, items2 = zip(*comb) if comb else ((), ())
                d1 = tuple(i for i in items1 if i)
                d2 = tuple(i for i in items2 if i)
                data.add((d1, d2))
        return DualMaySST2(data)

    def diff(self):
        """ return the boundary of the chain """
        data = set()
        for mon in self.data:
            for item1, item2 in itertools.combinations(mon, 2):
                k1, r1 = item1
                k2, r2 = item2
                i1, j1 = k1
                i2, j2 = k2
                if i1 + j1 == i2 or i2 + j2 == i1:
                    k = min(i1, i2), j1 + j2
                    mon_diff = dict(mon)
                    if k not in mon_diff or mon_diff[k] % 2 == 0:
                        mon_diff[k1] -= 1
                        mon_diff[k2] -= 1
                        if k in mon_diff:
                            mon_diff[k] += 1
                        else:
                            mon_diff[k] = 1
                        data ^= {tuple(sorted(item for item in mon_diff.items() if item[1]))}
        return DualMaySS(data)

    def inv_diff(self) -> bool:
        s, t, u = self.deg()
        my_map2 = linalg.LinearMapKernelMod2()
        my_map2.add_maps((r, r.diff()) for r in self.basis(s + 1, t, u))
        return my_map2.g(self)

    @classmethod
    def homology(cls, s, t, u):
        my_map1 = linalg.LinearMapKernelMod2()
        my_map2 = linalg.LinearMapKernelMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        print("kernel dim:", my_map1.kernel.dim())
        # for r in my_map1.kernel.basis(DualMaySS):
        #     print(r)
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s + 1, t, u))
        print("image: dim", my_map2.image().dim())
        # for r in my_map2.image().basis(DualMaySS):
        #     print(r)
        print("quotient:")
        result = list(my_map1.kernel.quotient(my_map2.image()).simplify().basis(DualMaySS))
        for r in result:
            print(r)
        return result

    def is_primitive(self):
        """ assert self.diff() == 0 """
        coprod = self.coprod()
        my_dict = {}
        s = self.deg()[0]
        for m1, m2 in coprod.data:
            if m1 and m2:
                tu1, tu2 = self.deg_mon(m1)[1:], self.deg_mon(m2)[1:]
                if (tu1, tu2) in my_dict:
                    my_dict[(tu1, tu2)] ^= {(m1, m2)}
                else:
                    my_dict[(tu1, tu2)] = {(m1, m2)}
        for tu1, tu2 in sorted(my_dict):
            t1, u1 = tu1
            t2, u2 = tu2
            if any((i, t1, u1) in self._homology and (s - i, t2, u2) in self._homology for i in range(s + 1)):
                boundaries = linalg.VectorSpaceMod2()
                for i in range(s + 2):
                    for m1 in self.basis_mons(i, t1, u1):
                        for m2 in self.basis_mons(s + 1 - i, t2, u2):
                            boundaries.add_v(DualMaySST2((m1, m2)).diff())
                if boundaries.res(my_dict[(tu1, tu2)]):
                    # print("${}$\\\\".format(DualMaySST2(my_dict[(tu1, tu2)])))
                    return False
        return True

    @classmethod
    def homology_basis(cls) -> Iterable["DualMaySS"]:
        return map(DualMaySS, itertools.chain.from_iterable(cls._homology.values()))

    @classmethod
    def search_primitives(cls):
        result = sorted((r for r in cls.homology_basis() if r.is_primitive()), key=lambda x: x.deg())
        for r in result:
            print(r)

    @staticmethod
    def lexicographic(mon) -> list:
        return sorted(mon, key=lambda g: (g[0][0] + g[0][1], -g[0][0]))


class DualMaySST2(BC.AlgebraT2Mod2):
    """ Tensor product of two DualSteenrod """
    type_c0 = DualMaySS
    type_c1 = DualMaySS

    def mul_mons(self, mon1: tuple, mon2: tuple):
        prod0 = self.type_c0.mul_mons(mon1[0], mon2[0])
        prod1 = self.type_c1.mul_mons(mon1[1], mon2[1])
        if type(prod0) is type(prod1) is tuple:
            return prod0, prod1
        else:
            return set()

    @classmethod
    def unit(cls):
        return cls(((), ()))

    def diff(self):
        data = set()
        for m1, m2 in self.data:
            for m in DualMaySS(m1).diff().data:
                data.add((m, m2))
            for m in DualMaySS(m2).diff().data:
                data.add((m1, m))
        return type(self)(data)

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (
            self.deg_mon(m), sorted(m[0].items()), sorted(m[1].items())), reverse=True)


def test():
    sig = Signature((1, 1, 0, -1, -1)).accumulate()
    s = 3
    basis_s = set()
    leading_terms = set()
    for m in DualMaySS.basis_sig_mon(sig):
        length = DualMaySS.deg_s_mon(m)
        if length == s:
            basis_s.add(m)
        elif length == s + 1:
            m = max(DualMaySS(m).diff().data, key=key_lex)
            if m is not None:
                leading_terms.add(m)
    for m in leading_terms:
        print(m)
    return basis_s - leading_terms


def key_lex(mon):
    return sorted(map(lambda g: (g[0][0] + g[0][1], -g[0][0]), mon))


if __name__ == "__main__":
    print("MaySS")

# 389
