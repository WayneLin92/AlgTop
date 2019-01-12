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
        for r in my_map1.kernel.get_basis(MaySS):
            print(r)
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s - 1, t, u))
        print("image: dim", my_map2.get_image().get_dim())
        for r in my_map2.get_image().get_basis(MaySS):
            print(r)
        print("quotient:")
        for r in my_map1.kernel.quotient(my_map2.get_image()).get_basis(MaySS):
            print(r)


class DualMaySS(BC.BaseExteriorMod2):
    """ This is the dual of the May spectral sequence """
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *key: int):
        assert len(key) == 3 and key[0] > 0 <= key[1]
        return cls(frozenset((ij2deg(key), 1 << k) for k in mymath.two_expansion(key[2])))

    @staticmethod
    def deg_gen(key: tuple) -> int:
        return key[0] * key[1]

    @staticmethod
    def str_gen(key: tuple) -> str:
        deg, n = key
        i, j = deg2ij(deg)
        return "\\gamma_{}(\\bar{{P}}_{}^{})".format(*map(mymath.tex_index, (n, i, j)))

    @classmethod
    def str_mon(cls, mon: frozenset):
        d = {}
        for deg, n in mon:
            if deg in d:
                d[deg] += n
            else:
                d[deg] = n
        result = "".join(map(cls.str_gen, d.items()))
        return result if result else "1"

    # methods -----------------
    @staticmethod
    def basis_mons(length, deg, may_filtr, ij_max=None):
        """ return monomials h_ij^k...h_index_max^e with given length, deg, may_filtr """
        cls = DualMaySS
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
                yield frozenset({(ij2deg(ij_max), deg)})
            return
        for e in range(min(length, deg // ij2deg(ij_max), may_filtr // ij_max[0]), -1, -1):
            # print(length, deg, may_filtr, "{}^{}".format(ij_max, e))
            index_next = (ij_max[0] - 1, ij_max[1] + 1) if ij_max[0] > 1 else (sum(ij_max) - 1, 0)
            for mon in cls.basis_mons(length - e, deg - e * ij2deg(ij_max),
                                      may_filtr - e * ij_max[0], index_next):
                if e > 0:
                    yield mon ^ frozenset({(ij2deg(ij_max), e)})
                else:
                    yield mon

    @classmethod
    def basis(cls, length, deg, may_filtr):
        return (cls(m) for m in DualMaySS.basis_mons(length, deg, may_filtr))

    def diff(self):
        pass


# functions
def ij2deg(key) -> int:
    return (1 << key[0]) - 1 << key[1]


def deg2ij(n: int):
    i = bin(n).count('1')
    j = n.bit_length() - i
    return i, j
