import itertools
from algebras import BaseClasses as BC, linalg, mymath


class MaySS(BC.BasePolyMod2):
    """ This is for the MSS starting from E_1 """
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *key: int):
        assert len(key) == 2 and key[0] > 0 <= key[1]
        return cls(((ij2deg(key), 1),))

    @staticmethod
    def deg_gen(n: int) -> int:
        return n

    @staticmethod
    def str_gen(n: int) -> str:
        return "h_{{{}，{}}}".format(*deg2ij(n))

    # methods -----------------
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
                yield ((ij2deg(ij_max), deg),)
            return
        for e in range(min(length, deg // ij2deg(ij_max), may_filtr // ij_max[0]), -1, -1):
            # print(length, deg, may_filtr, "{}^{}".format(ij_max, e))
            index_next = (ij_max[0] - 1, ij_max[1] + 1) if ij_max[0] > 1 else (sum(ij_max) - 1, 0)
            for mon in cls.basis_mons(length - e, deg - e * ij2deg(ij_max),
                                      may_filtr - e * ij_max[0], index_next):
                if e > 0:
                    yield mon + ((ij2deg(ij_max), e),)
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
                    i, j = deg2ij(g)
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
        print("kernel dim:", my_map1.kernel.get_dim())
        # for r in my_map1.kernel.get_basis(MaySS):
        #     print(r)
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s - 1, t, u))
        print("image: dim", my_map2.get_image().get_dim())
        # for r in my_map2.get_image().get_basis(MaySS):
        #     print(r)
        print("quotient:")
        for r in my_map1.kernel.quotient(my_map2.get_image()).get_basis(MaySS):
            print(r)


class DualMaySS(BC.BaseExteriorMod2):
    """ This is the dual of the May spectral sequence """
    _maps = None
    _homology = None

    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *key: int):
        assert len(key) == 3 and key[0] > 0 <= key[1]
        return cls(frozenset((ij2deg(key), 1 << k) for k in mymath.two_expansion(key[2])))

    @staticmethod
    def deg_gen(gen: tuple) -> mymath.Deg:
        deg, r = gen
        i, j = deg2ij(deg)
        return mymath.Deg((r, deg * r, i * r))

    @staticmethod
    def str_gen(gen: tuple) -> str:
        deg, n = gen
        i, j = deg2ij(deg)
        return "\\gamma_{}(\\bar{{P}}_{}^{})".format(*map(mymath.tex_index, (n, i, j)))

    @classmethod
    def str_mon(cls, mon: frozenset):
        result = "".join(map(cls.str_gen, cls.comb_gens(mon)))
        return result if result else "1"

    @classmethod
    def deg_mon(cls, mon: frozenset):
        return super().deg_mon(mon) if mon else mymath.Deg((0, 0, 0))

    # methods -----------------
    @staticmethod
    def basis_mons(length, deg, may_filtr, ij_max=None):
        """ return monomials h_ij^k...h_index_max^e with given length, deg, may_filtr """
        if length == 0 or deg == 0 or may_filtr == 0:
            if length == deg == may_filtr:
                yield frozenset()
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
                yield frozenset((ij2deg(ij_max), 1 << k) for k in mymath.two_expansion(deg))
            return
        for e in range(min(length, deg // ij2deg(ij_max), may_filtr // ij_max[0]), -1, -1):
            index_next = (ij_max[0] - 1, ij_max[1] + 1) if ij_max[0] > 1 else (sum(ij_max) - 1, 0)
            for mon in DualMaySS.basis_mons(length - e, deg - e * ij2deg(ij_max),
                                            may_filtr - e * ij_max[0], index_next):
                if e > 0:
                    yield mon | frozenset((ij2deg(ij_max), 1 << k) for k in mymath.two_expansion(e))
                else:
                    yield mon

    @classmethod
    def basis(cls, length, deg, may_filtr):
        return (cls(m) for m in DualMaySS.basis_mons(length, deg, may_filtr))

    @staticmethod
    def comb_gens(mon: frozenset):
        """ return an iterator of (gen, r) with r's combined """
        mon_dict = {}
        for deg, n in mon:
            if deg in mon_dict:
                mon_dict[deg] += n
            else:
                mon_dict[deg] = n
        return mon_dict.items()

    @staticmethod
    def coprod_gen(gen):
        deg, r = gen
        data = {(frozenset((deg, 1 << e) for e in mymath.two_expansion(k)),
                 frozenset((deg, 1 << e) for e in mymath.two_expansion(r - k)))
                for k in range(r + 1)}
        return DualMaySST2(data)

    def coprod(self):
        result = DualMaySST2.zero()
        for m in self.data:
            product = result.unit()
            for gen in self.comb_gens(m):
                product = product * self.coprod_gen(gen)
            result += product
        return result

    def diff(self):
        """ return the boundary of the chain """
        result = self.zero()
        for mon in self.data:
            mon_min_gamma = {}
            for deg, n in mon:
                if deg in mon_min_gamma:
                    if mon_min_gamma[deg] > n:
                        mon_min_gamma[deg] = n
                else:
                    mon_min_gamma[deg] = n
            for gs, gt in itertools.combinations(mon_min_gamma.items(), 2):
                deg_s, r_s = gs
                deg_t, r_t = gt
                i_s, j_s = deg2ij(deg_s)
                i_t, j_t = deg2ij(deg_t)
                if j_s == i_t + j_t or j_t == i_s + j_s:
                    factor1 = type(self).gen(i_s + i_t, min(j_s, j_t), 1)
                    gs_minus = set((deg_s, 1 << k) for k in range(r_s.bit_length() - 1))
                    gt_minus = set((deg_t, 1 << k) for k in range(r_t.bit_length() - 1))
                    factor2 = type(self)(mon - {gs, gt} | gs_minus | gt_minus)
                    result += factor1 * factor2
        return result

    @classmethod
    def homology(cls, s, t, u):
        my_map1 = linalg.LinearMapKernelMod2()
        my_map2 = linalg.LinearMapKernelMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        print("kernel dim:", my_map1.kernel.get_dim())
        for r in my_map1.kernel.get_basis(DualMaySS):
            print(r)
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s + 1, t, u))
        print("image: dim", my_map2.get_image().get_dim())
        for r in my_map2.get_image().get_basis(DualMaySS):
            print(r)
        print("quotient:")
        for r in my_map1.kernel.quotient(my_map2.get_image()).get_basis(DualMaySS):
            print(r)

    @classmethod
    def load(cls, s_max, t_max, u_max):
        # TODO: load from file
        cls._maps = {}
        cls._homology = {}
        for s in range(s_max + 2):
            for t in range(s, t_max + 1):
                for u in range(s, u_max + 1):
                    cls._maps[(s, t, u)] = linalg.LinearMapKernelMod2()
                    cls._maps[(s, t, u)].add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        for s in range(s_max + 1):
            for t in range(s, t_max + 1):
                for u in range(s, u_max + 1):
                    cycles = cls._maps[(s, t, u)].kernel
                    if (s + 1, t, u) in cls._maps:
                        boundaries = cls._maps[(s + 1, t, u)].get_image()
                        cls._homology[(s, t, u)] = list(cycles.quotient(boundaries).simplify().get_basis(DualMaySS))
                    else:
                        cls._homology[(s, t, u)] = list(cycles.simplify().get_basis(DualMaySS))
                    if not cls._homology[(s, t, u)]:
                        del cls._homology[(s, t, u)]

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
        for tu1, tu2 in my_dict:
            t1, u1 = tu1
            t2, u2 = tu2
            if any((i, t1, u1) in self._homology and (s - i, t2, u2) in self._homology for i in range(s + 1)):
                print(tu1, tu2)
                boundaries = linalg.VectorSpaceMod2()
                for i in range(s + 2):
                    for m1 in self.basis_mons(i, t1, u1):
                        for m2 in self.basis_mons(s + 1 - i, t2, u2):
                            boundaries.add_v(DualMaySST2((m1, m2)).diff())
                if boundaries.res(my_dict[(tu1, tu2)]):
                    return False
        return True

    @classmethod
    def search_primitives(cls):
        for r in cls.homology_basis():
            if r.is_primitive():
                print("${}$\\\\".format(r))

    @classmethod
    def homology_basis(cls):
        return itertools.chain.from_iterable(cls._homology.values())


class DualMaySST2(BC.GradedRingT2Mod2):
    """ Tensor product of two DualSteenrod """
    type_c0 = DualMaySS
    type_c1 = DualMaySS

    def mul_mons(self, mon1: tuple, mon2: tuple):
        prod0 = self.type_c0.mul_mons(mon1[0], mon2[0])
        prod1 = self.type_c1.mul_mons(mon1[1], mon2[1])
        if type(prod0) is type(prod1) is frozenset:
            return prod0, prod1
        else:
            return set()

    @classmethod
    def unit(cls):
        return cls((frozenset(), frozenset()))

    def diff(self):
        data = set()
        for m1, m2 in self.data:
            for m in DualMaySS(m1).diff().data:
                data.add((m, m2))
            for m in DualMaySS(m2).diff().data:
                data.add((m1, m))
        return type(self)(data)


# functions
def ij2deg(key: tuple) -> int:
    return (1 << key[0]) - 1 << key[1]


def deg2ij(n: int) -> tuple:
    i = bin(n).count('1')
    j = n.bit_length() - i
    return i, j


if __name__ == "__main__":
    DualMaySS.load(5, 30, 20)
    L = list(itertools.chain.from_iterable(DualMaySS._homology.values()))
    for rr in L:
        if len(rr.data) > 1:
            print("${}: {}$\\\\".format(rr.deg(), rr))
