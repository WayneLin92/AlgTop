import itertools
import heapq
from algebras import BaseAlgebras as BA, mymath


class AlgB(BA.BasePolyMod2):
    _rels = {}
    _index_max = 0

    # ----- AlgebraMod2 -------------
    @classmethod
    def gen(cls, *key: int) -> "AlgB":
        """Return (B^i_j)^r."""
        assert key[0] < key[1] <= cls._index_max
        if len(key) == 2:
            return cls(((key, 1),))
        else:
            i, j, r = key
            return cls((((i, j), r),))

    @staticmethod
    def deg_gen(k) -> mymath.Deg:
        i, j = k
        return mymath.Deg((1, 2 ** j - 2 ** i, j - i)) * 2

    @staticmethod
    def str_gen(k) -> str:
        return "B^{}_{}".format(*map(mymath.tex_index, k))

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

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=self.key_mon)

    # methods -----------------
    @classmethod
    def set_index_max(cls, i_max):
        cls._index_max = i_max
        B = cls.gen
        for j in range(2, i_max + 1):
            for i in range(i_max - j + 1):
                s, t = i, i + j
                rel = sum((B(s, k) * B(k, t) for k in range(s + 1, t)), cls.zero())
                cls._add_rel(rel)

    @staticmethod
    def key_mon(mon):
        n_cross = sum(g1[1] * g2[1] for g1, g2 in itertools.combinations(mon, 2)
                      if g1[0][0] < g2[0][0] < g1[0][1] < g2[0][1])
        return -n_cross, [(k, -r) for k, r in mon]

    @classmethod
    def leading(cls, data):
        """Return the leading monomial of data."""
        return min(data, key=cls.key_mon)

    @staticmethod
    def gcd_nonzero(mon1, mon2):
        """Return gcd(mon1, mon2) > 0."""
        keys1 = {k for k, r in mon1}
        return any(k in keys1 for k, r in mon2)

    @staticmethod
    def divisible(mon1, mon2):
        """Return mon1 | mon2."""
        d2 = dict(mon2)
        return all(k in d2 and r <= d2[k] for k, r in mon1)

    @staticmethod
    def div_mons(mon1, mon2):
        """Return mon1 / mon2."""
        d2 = dict(mon2)
        return tuple((k, r - d2[k] if k in d2 else r) for k, r in mon1 if k not in d2 or r > d2[k])

    @staticmethod
    def lcm_mons(mon1, mon2):
        """Return lcm(mon1, mon2)."""
        result = dict(mon1)
        for k, r in mon2:
            result[k] = max(result[k], r) if k in result else r
        return tuple(sorted(result.items()))

    @staticmethod
    def div_mod(mon1, mon2):
        """Return q := exp1 // exp2, r := mon1 - mon2 ** q."""
        d2 = dict(mon2)
        q = min(r // d2[k] for k, r in mon1 if k in d2)
        if q:
            mon2_q = tuple((k, r * q) for k, r in mon2)
            return q, AlgB.div_mons(mon1, mon2_q)
        else:
            return q, ()

    @classmethod
    def _add_rel(cls, rel):
        """Add a relation."""
        print("  rel:", cls(rel) if type(rel) is set else rel)
        if not rel:
            return
        if type(rel) is not set:
            hq = [(rel.deg()[1], rel.data)]
        else:
            hq = [(cls(rel).deg()[1], rel)]
        while hq:
            deg, r = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                m = cls.leading(r)
                print("    leading:", cls(m), "deg:", deg)
                redundant_leading_terms = []
                for m1, v1 in cls._rels.items():
                    if cls.gcd_nonzero(m, m1):  # gcd > 0
                        if cls.divisible(m, m1):  # m | m1
                            redundant_leading_terms.append(m1)
                            heapq.heappush(hq, (cls.deg_mon(m1)[1], v1 | {m1}))
                        else:
                            lcm = cls.lcm_mons(m, m1)
                            dif = cls.div_mons(lcm, m)
                            dif1 = cls.div_mons(lcm, m1)
                            new_rel = {cls.mul_mons(_m, dif) for _m in r}
                            v1dif1 = {cls.mul_mons(_m, dif1) for _m in v1}
                            new_rel -= {lcm}
                            new_rel ^= v1dif1
                            heapq.heappush(hq, (cls.deg_mon(lcm)[1], new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls._rels[m_redundant]
                cls._rels[m] = r - {m}

    @classmethod
    def simplify_rels(cls):
        """Simplify cls._rels.
        Should be called after the completion of the construction of the algebra.
        """
        for m in cls._rels:
            cls._rels[m] = cls.simplify_data(cls._rels[m])

    @classmethod
    def simplify_data(cls, data: set):
        """Simplify the data by relations."""
        s = data.copy()
        result = set()
        while s:
            mon = cls.leading(s)
            s.remove(mon)
            for m in cls._rels:
                if cls.divisible(m, mon):
                    q, r = cls.div_mod(mon, m)
                    s ^= {cls.mul_mons(r, tuple((k, r * q) for k, r in m1))
                          for m1 in cls._rels[m]}
                    break
            else:
                if mon in result:  # ##############
                    print("mon =", cls(mon))
                    print("s =", cls(s))
                    assert False
                result.add(mon)
        return result

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    @classmethod
    def print_rels(cls):
        for m in sorted(cls._rels):
            if any(k[0] == 0 for k, r in m):
                print(f"${cls(m)} = {cls(cls._rels[m])}$\\\\")


if __name__ == "__main__":
    AlgB.set_index_max(7)
    AlgB.print_rels()
