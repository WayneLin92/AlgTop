import itertools
import operator
import pickle
from algebras import BaseAlgebras as BC, linalg, mymath


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
        return Signature(map(operator.mul, self, itertools.repeat(other)))

    def __bool__(self):
        """Return if any element in self is nonzero."""
        return any(self)

    def simplify(self) -> "Signature":
        """Remove trailing zeroes."""
        i = len(self)
        while i > 0 and self[i-1] == 0:
            i -= 1
        return Signature(self[:i])

    def deg(self):
        """Return the degree of the monomial with self as the signature."""
        return sum(fn << n for n, fn in enumerate(self))

    def span(self):
        """Return [left, r) on which self is non-zero."""
        right = len(self)
        while right > 0 and self[right-1] == 0:
            right -= 1
        left = 0
        while left < right and self[left] == 0:
            left += 1
        return left, right

    def accumulate(self):
        """Return the sequence of the partial sum."""
        return Signature(itertools.accumulate(self[:-1]))

    def diff(self):
        """Return the sequence of self[i] - self[i-1]."""
        return Signature((self[0],) + tuple(self[i] - self[i-1] for i in range(1, len(self))) + (-self[-1],))


# noinspection PyUnresolvedReferences
class MayDGA:
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
                for m in MayDGA._fixed_i_sig_mon(sig - ((0,) * i + (r,) * j), i, j - 1):
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
        for m in MayDGA._fixed_i_sig_mon(sig, i):
            sig_m = MayDGA.sig_mon(m)
            for m1 in MayDGA.basis_sig_mon(sig - sig_m):
                yield m + m1

    @classmethod
    def basis_sig(cls, sig: Signature):
        return map(cls, cls.basis_sig_mon(sig))

    @classmethod
    def homology_signature(cls, sig: Signature):
        s_min = sum(i for i in sig.diff() if i > 0)
        s_max = sum(sig)
        lin_maps = {}
        for s in range(s_min - 1, s_max + 2):
            lin_maps[s] = linalg.LinearMapKMod2()
        for x in cls.basis_sig(sig):
            s = x.deg()[0]
            lin_maps[s].add_maps(((x, x.diff()),))
        homology = {}
        for s in range(s_min, s_max + 1):
            print(f"{s}:")
            s1 = s - 1 if cls is MaySS else s + 1
            homology[s] = list((lin_maps[s].kernel / lin_maps[s1].image).simplify().basis(cls))
            for x in homology[s]:
                print(x)
        return homology


class MaySS(MayDGA, BC.BasePolyMod2):
    """This is for the May spectral sequence starting from E_1."""
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *key: int):
        """Return (R^i_j)^r."""
        if len(key) == 2:
            i, j = key
            return cls((((i, j), 1),))
        else:
            i, j, r = key
            return cls((((i, j), r),))

    @staticmethod
    def deg_gen(k) -> mymath.Deg:
        i, j = k
        return mymath.Deg((1, (1 << j) - 1 << i, j))

    @staticmethod
    def str_gen(k):
        return "K^{}_{}".format(*map(mymath.tex_index, (k[0], k[0] + k[1])))

    @classmethod
    def str_mon(cls, mon):
        result = ""
        for gen, exp in mon:
            if exp == 1:
                result += cls.str_gen(gen)
            else:
                result += f"({cls.str_gen(gen)})^{mymath.tex_index(exp)}"
        if result == "":
            result = "1"
        return result

    # methods -----------------
    @classmethod
    def h(cls, i: int, S: tuple):
        """Return h_i(S)."""
        k = len(S)
        seq = set(range(i, i + 2 * k + 2))
        S = {i + s for s in S} | {i}
        T = seq - S
        assert len(S) + len(T) == 2 * k + 2
        S, T = sorted(S), sorted(T)
        data = set()
        for T1 in itertools.permutations(T):
            J = tuple(map(operator.sub, T1, S))
            if all(j > 0 for j in J):
                m = tuple(zip(zip(S, J), itertools.repeat(1)))
                data.add(m)
        return cls(data)

    @classmethod
    def Phi(cls, S, T):
        """Return H^S_T."""
        assert len(S) == len(T) == len(set(S)) == len(set(T))
        S, T = sorted(S), sorted(T)
        data = set()
        for T1 in itertools.permutations(T):
            J = tuple(map(operator.sub, T1, S))
            if all(j > 0 for j in J):
                m = tuple(zip(zip(S, J), itertools.repeat(1)))
                data.add(m)
        return cls(data)

    @staticmethod
    def deg_t_gen(k: tuple) -> int:
        return (1 << k[1]) - 1 << k[0]

    @staticmethod
    def basis_mons(deg_s, deg_t, deg_u, ij_max=None):
        """Return monomials (R^i_j)^k...(R ij_max)^e with given length, deg, may_filtr."""
        if deg_s == 0 or deg_t == 0 or deg_u == 0:
            if deg_s == deg_t == deg_u:
                yield ()
            return
        if ij_max is None:
            bound = deg_t - deg_s + 1
            if bound <= 0:
                return
            i_plus_j = bound.bit_length()
            i = 0
            while (1 << i_plus_j) - (1 << i) > bound:
                i += 1
            ij_max = (i, i_plus_j - i)
        if ij_max == (0, 1):
            if deg_s == deg_t == deg_u:
                yield (((0, 1), deg_t),)
            return
        i, j = ij_max
        for e in range(min(deg_s, deg_t // MaySS.deg_t_gen(ij_max), deg_u // j), -1, -1):
            ij_next = (i + 1, j - 1) if j > 1 else (0, i + j - 1)
            for mon in MaySS.basis_mons(deg_s - e, deg_t - e * MaySS.deg_t_gen(ij_max),
                                        deg_u - e * j, ij_next):
                if e > 0:
                    yield tuple(sorted(mon + ((ij_max, e),)))
                else:
                    yield mon

    @classmethod
    def basis(cls, deg_s, deg_t, deg_u):
        return map(cls, MaySS.basis_mons(deg_s, deg_t, deg_u))

    def diff(self):
        """Return the coboundary of the cochain."""
        result = self.zero()
        for m in self.data:
            for ind in range(len(m)):
                k, e = m[ind]
                if e % 2:
                    i, j = k
                    dk_data = set((((i, j1), 1), ((i + j1, j - j1), 1)) for j1 in range(1, j))
                    dk = MaySS(dk_data)
                    if e - 1:
                        result += MaySS(m[:ind] + m[ind + 1:]) * MaySS(((k, e - 1),)) * dk
                    else:
                        result += MaySS(m[:ind] + m[ind + 1:]) * dk
        return result

    def inv_diff(self) -> bool:
        s, t, u = self.deg()
        my_map2 = linalg.LinearMapKMod2()
        my_map2.add_maps((r, r.diff()) for r in self.basis(s - 1, t, u))
        return my_map2.g(self)

    @classmethod
    def homology(cls, s, t, u):
        my_map1 = linalg.LinearMapKMod2()
        my_map2 = linalg.LinearMapKMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        print("kernel dim:", my_map1.kernel.dim)
        s1 = s - 1 if cls is MaySS else s + 1
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s1, t, u))
        print("image: dim", my_map2.image.dim)
        print("quotient:")
        for r in (my_map1.kernel / my_map2.image).basis(cls):
            print(r)

    def print_tex_graph(self):
        print_tex_graph(self.data, sep="\\hspace{5pt}+\\hspace{5pt}")


class DualMaySS(MayDGA, BC.AlgebraMod2):
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
                        cls._maps[(s, t, u)] = linalg.LinearMapKMod2()
                        cls._maps[(s, t, u)].add_maps((r, r.diff()) for r in cls.basis(s, t, u))
                        if s <= s_max:
                            cycles = cls._maps[(s, t, u)].kernel
                            if (s + 1, t, u) in cls._maps:
                                boundaries = cls._maps[(s + 1, t, u)].image
                                cls._homology[(s, t, u)] = list((cycles / boundaries).simplify().basis(set))
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
        return cls((((i, j), r),))

    @staticmethod
    def str_item(item: tuple) -> str:
        ij, r = item
        return "\\gamma_{}(\\bar{{P}}^{}_{})".format(*map(mymath.tex_index, (r, *ij)))

    @classmethod
    def str_mon(cls, mon: tuple):
        mon = cls.lexicographic(mon)
        result = "".join(map(cls.str_item, mon))
        return result if result else "1"

    @staticmethod
    def deg_item(item: tuple):
        ij, r = item
        return mymath.Deg((r, ((1 << ij[1]) - 1 << ij[0]) * r, ij[1] * r))

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(cls.deg_item, mon), mymath.Deg((0, 0, 0)))

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (self.deg_mon(m), m))

    # methods -----------------
    deg_t_k = MaySS.deg_t_gen
    basis_mons = MaySS.basis_mons
    basis = MaySS.basis
    homology = MaySS.homology

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
    def homology_basis(cls):
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


def key_lex(mon):
    return sorted(map(lambda g: (-g[0][0] - g[0][1], g[0][0]), mon))


def key_lex_reverse(mon):
    return sorted(map(lambda g: (g[0][0] + g[0][1], -g[0][0]), mon), reverse=True)


def tex_graph_mon(mon):
    sig = DualMaySS.sig_mon(mon)
    left, right = sig.span()
    sep = 6
    result = f"\\xymatrix@M=0pt@C={sep}pt{{\n"
    arrows = [[] for _ in range(right-left)]
    for k, r in mon:
        # noinspection PyTypeChecker
        arrows[k[0] - left].append((k[1], r))
    for i in range(right - left):
        result += "\\bullet"
        for j, right in arrows[i]:
            result += f" \\ar@/^{sep*j}pt/@{{-}}[{'r'*j}]"
            if right > 1:
                result += f"|-{right}"
        result += " & "
    result += "\\bullet\n}"
    return result


def print_tex_graph(iterable, *, row=5, sep=',\\hspace{5pt}'):
    i = 0
    for i, m in enumerate(iterable):
        if i % row == 0:
            print("$$")
        print(f"{tex_graph_mon(m)}{sep}")
        if i % row == row - 1:
            print("$$\n")
    if i % row != row - 1:
        print("$$\n")


def Phi(S, T=None):
    if type(S) is str:
        nums = [int(c) for c in S if c in "0123456789"]
        assert len(nums) % 2 == 0
        half = len(nums) // 2
        return MaySS.Phi(nums[:half], nums[half:])
    elif type(S) is int:
        return MaySS.Phi((S,), [T])
    else:
        return MaySS.Phi(S, T)


def test():
    pass


if __name__ == "__main__":
    pass

# 389, 551, 596, 536, 542
