import itertools
import operator

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
        """Return monomials of gamma_r P^i_j for fixed i"""
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
            s1 = s - 1 if cls is MayE1 else s + 1
            homology[s] = list((lin_maps[s].kernel / lin_maps[s1].image).simplify().basis(cls))
            for x in homology[s]:
                print(x)
        return homology


class MayE1(MayDGA, BC.BasePolyMod2):
    """E_1 page of the May spectral sequence."""
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *k: int):
        """Return R_{ij}."""
        i, j = k
        return cls((((i, j), 1),))

    @staticmethod
    def deg_gen(k) -> mymath.Deg:
        i, j = k
        return mymath.Deg((1, (1 << j) - (1 << i), j - i))

    @classmethod
    def deg_mon(cls, mon: frozenset) -> mymath.Deg:
        return sum((cls.deg_gen(gen) * exp for gen, exp in mon), mymath.Deg((0, 0, 0)))

    @staticmethod
    def str_gen(k):
        return f"R_{{{k[0]}{k[1]}}}"

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
            if all(t - s > 0 for s, t in zip(S, T1)):
                m = tuple(zip(zip(S, T1), itertools.repeat(1)))
                data.add(m)
        return cls(data)

    @staticmethod
    def deg_t_gen(k: tuple) -> int:
        return (1 << k[1]) - (1 << k[0])

    @staticmethod
    def basis_mons(deg_s, deg_t, deg_u, ij_max=None):
        """Return monomials (R_{i_1j_1})^k...(R i_mj_m)^e with given length, deg, may_filtr."""
        if deg_s == 0 or deg_t == 0 or deg_u == 0:
            if deg_s == deg_t == deg_u:
                yield ()
            return
        if ij_max is None:
            bound = deg_t - deg_s + 1
            if bound <= 0:
                return
            j = bound.bit_length()
            i = 0
            while (1 << j) - (1 << i) > bound:
                i += 1
            ij_max = (i, j)
        if ij_max == (0, 1):
            if deg_s == deg_t == deg_u:
                yield (((0, 1), deg_t),)
            return
        i, j = ij_max
        for e in range(min(deg_s, deg_t // MayE1.deg_t_gen(ij_max), deg_u // (j - i)), -1, -1):
            ij_next = (i + 1, j) if j - i > 1 else (0, j - 1)
            for mon in MayE1.basis_mons(deg_s - e, deg_t - e * MayE1.deg_t_gen(ij_max),
                                        deg_u - e * (j - i), ij_next):
                if e > 0:
                    yield tuple(sorted(mon + ((ij_max, e),)))
                else:
                    yield mon

    @classmethod
    def basis(cls, deg_s, deg_t, deg_u):
        return map(cls, MayE1.basis_mons(deg_s, deg_t, deg_u))

    def diff(self):
        """Return the differential."""
        result = self.zero()
        for m in self.data:
            for ind in range(len(m)):
                k, e = m[ind]
                if e % 2:
                    i, j = k
                    dk_data = set((((i, _k), 1), ((_k, j), 1)) for _k in range(i + 1, j))
                    dk = MayE1(dk_data)
                    if e - 1:
                        result += MayE1(m[:ind] + m[ind + 1:]) * MayE1(((k, e - 1),)) * dk
                    else:
                        result += MayE1(m[:ind] + m[ind + 1:]) * dk
        return result

    @classmethod
    def homology(cls, s, t, u):
        my_map1 = linalg.LinearMapKMod2()
        my_map2 = linalg.LinearMapKMod2()
        my_map1.add_maps((r, r.diff()) for r in cls.basis(s, t, u))
        print("kernel dim:", my_map1.kernel.dim)
        s1 = s - 1 if cls is MayE1 else s + 1
        my_map2.add_maps((r, r.diff()) for r in cls.basis(s1, t, u))
        print("image dim:", my_map2.image.dim)
        print("quotient:")
        for r in (my_map1.kernel / my_map2.image).basis(cls):
            print(r)

    def _repr_latex_(self):
        return f"${self}$"


class TorMayE1(BC.BaseExteriorMod2):
    """Class for Tor_{MayE1}(F_2, F_2)."""
    # ----- BasePolyMod2 ----------------
    @classmethod
    def gen(cls, *k: int):
        """Return hat R_{ij}."""
        i, j = k
        return cls(frozenset({(i, j)}))

    @staticmethod
    def deg_gen(k) -> mymath.Deg:
        i, j = k
        return mymath.Deg((1, (1 << j) - (1 << i), j - i))

    @classmethod
    def deg_mon(cls, mon: frozenset) -> mymath.Deg:
        return sum(map(cls.deg_gen, mon), mymath.Deg((0, 0, 0)))

    @staticmethod
    def str_gen(k) -> str:
        return f"\\hat R_{{{k[0]}{k[1]}}}"

    @staticmethod
    def basis_mons(n_max):
        """Return monomials (R_{i_1j_1})^k...(R i_mj_m)^e with given length, deg, may_filtr."""
        generators = []
        for i in range(n_max):
            for j in range(i + 1, n_max + 1):
                generators.append((i, j))
        for r in range(len(generators)+1):
            for m in itertools.combinations(generators, r):
                yield frozenset(m)


class MayE1Res(BC.AlgebraT2Mod2):
    """Resolution of MayE1."""
    type_c0 = TorMayE1
    type_c1 = MayE1

    def diff(self):
        """Return the differential."""
        result = self.zero()
        for m0, m1 in self.data:
            dm1 = self.type_c1(m1).diff()
            result += self.tensor(self.type_c0(m0), dm1)
            for i, j in m0:
                dk_data = {(frozenset({(i, k)}), (((k, j), 1),)) for k in range(i + 1, j)}
                dk_data |= {(frozenset(), (((i, j), 1),))}
                dk = MayE1Res(dk_data)
                result += MayE1Res((m0 - {(i, j)}, m1)) * dk
        return result

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (
            self.deg_mon(m), sorted(m[0]), m[1]), reverse=True)


def test():
    pass


if __name__ == "__main__":
    pass

# 389, 551, 596, 536, 542