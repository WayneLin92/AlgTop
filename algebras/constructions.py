"""Classes that can be constructed by other classes."""
# todo: create class AugModuleMod2
# Todo: construct AugAlgMod2 from other algebras
# todo: use {generator} instead of set(generator)
import copy
import itertools
import operator
import heapq
from typing import Union, Set, Tuple, List, Dict, Type
from algebras import BaseAlgebras as BA, linalg, mymath


class AugAlgMod2(BA.AlgebraMod2):
    """A factory for commutative augmented graded algebras over F_2.

    AugAlgMod2.new_alg() creates a new commutative augmented graded algebra with
    its own generators and relations.
    """

    _gen_names = None  # type: List[str]
    _gen_degs = None  # type: list
    _unit_deg = None
    _rels = None  # type: Dict[tuple, set]
    _auto_simplify = None  # type: bool
    _name_index = 0

    @staticmethod
    def new_alg(unit_deg=None) -> "Type[AugAlgMod2]":
        """Return a dynamically created subclass of AugAlgMod2."""
        cls = AugAlgMod2
        class_name = f"AugAlgMod2_{cls._name_index}"
        cls._name_index += 1
        dct = {'_gen_names': [], '_gen_degs': [], '_unit_deg': unit_deg or 0, '_rels': {}, '_auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    # ----- AlgebraMod2 -------------
    @classmethod
    def mul_mons(cls, mon1: tuple, mon2: tuple):
        m = mymath.add_tuple(mon1, mon2)
        return cls.simplify_data({m}) if cls._auto_simplify else m

    @classmethod
    def str_mon(cls, mon: tuple):
        if mon:
            return "".join(f"{s}{mymath.tex_exponent(e)}" for e, s in zip(mon, cls._gen_names) if e)
        else:
            return "1"

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(operator.mul, mon, cls._gen_degs), cls._unit_deg)

    def deg(self):
        """self should always be homogeneous."""
        for m in self.data:
            return self.deg_mon(m)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, k: str, deg):
        """Add a new generator and return it."""
        cls._gen_names.append(k)
        cls._gen_degs.append(deg)
        m = (0,) * (len(cls._gen_names) - 1) + (1,)
        return cls(m).simplify()

    @classmethod
    def add_gens(cls, names, degs):
        """Add generators."""
        cls._gen_names += names
        cls._gen_degs += degs

    @classmethod
    def add_rel(cls, rel):
        """Add a relation."""
        if not rel:
            return
        if type(rel) is not set:
            hq = [(rel.deg(), rel.data)]
        else:
            hq = [(cls(rel).deg(), rel)]
        while hq:
            deg, r = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                m = min(r)
                redundant_leading_terms = []
                for m1, v1 in cls._rels.items():
                    if any(map(min, m, m1)):  # gcd > 0
                        if mymath.le_tuple(m, m1):
                            redundant_leading_terms.append(m1)
                            heapq.heappush(hq, (cls.deg_mon(m1), v1 | {m1}))
                        else:
                            lcm = mymath.max_tuple(m, m1)
                            dif = mymath.sub_tuple(lcm, m)
                            dif1 = mymath.sub_tuple(lcm, m1)
                            new_rel = {mymath.add_tuple(_m, dif) for _m in r}
                            v1dif1 = {mymath.add_tuple(_m, dif1) for _m in v1}
                            new_rel -= {lcm}
                            new_rel ^= v1dif1
                            heapq.heappush(hq, (cls.deg_mon(lcm), new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls._rels[m_redundant]
                cls._rels[m] = r - {m}
                # print(cls(m))

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
        s = list(data)
        heapq.heapify(s)
        result = set()

        leading_masks = tuple({i for i, e in enumerate(m) if e} for m in cls._rels)
        while s:
            mon = heapq.heappop(s)
            while s and mon == s[0]:
                heapq.heappop(s)
                mon = heapq.heappop(s) if s else None
            if mon is None:
                break
            mask_mon = {i for i, e in enumerate(mon) if e}
            for m, mask_m in zip(cls._rels, leading_masks):
                if mask_m <= mask_mon and mymath.le_tuple(m, mon):
                    q, r = mymath.div_mod_tuple(mon, m)
                    s += (mymath.add_tuple(r, tuple(map(operator.mul, m1, itertools.repeat(q))))
                          for m1 in cls._rels[m])
                    heapq.heapify(s)
                    break
            else:
                result.add(mon)
        return result

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    @staticmethod
    def _is_square(mon):
        return all(i % 2 == 0 for i in mon)

    @staticmethod
    def _first_non_square(data):
        for m in sorted(data):
            if any(i & 1 for i in m):
                return m

    @classmethod
    def _sqrt_zero(cls, lead: tuple, square: set, non_square: set, rels_lead_non_square: dict, deg_max_rel):
        """Generate an x^2(=0) where x is nonzero."""
        if not non_square:
            return square ^ {lead}
        mon = min(non_square)
        lcms1, lcms2 = {}, {}
        for m in cls._rels:
            if any(map(min, mon, m)):
                lcm = mymath.max_tuple(mon, m)
                dif = mymath.sub_tuple(lcm, mon)  # ############
                lcm = tuple(i + (j & 1) for i, j in itertools.zip_longest(lcm, dif, fillvalue=0))
                lcms1[lcm] = m
        for m in rels_lead_non_square:
            if any(map(min, mon, m)):
                lcm = mymath.max_tuple(mon, m)
                dif = mymath.sub_tuple(lcm, mon)  # ############
                lcm = tuple(i + (j & 1) for i, j in itertools.zip_longest(lcm, dif, fillvalue=0))
                dif1 = mymath.sub_tuple(lcm, m)
                if cls._is_square(dif1):
                    lcms2[lcm] = m

        lcms = [(m, 1) for m in lcms1]
        lcms += [(m, 2) for m in lcms2]
        lcms.sort(key=lambda _x: (cls.deg_mon(_x[0]), _x[1]))
        for lcm, index in lcms:
            dif = mymath.sub_tuple(lcm, mon)
            lead1 = mymath.add_tuple(lead, dif)
            lead1_half = tuple(i // 2 for i in lead1)
            if cls.is_admissible(lead1_half) and cls.deg_mon(lead1_half) < deg_max_rel:
                print(" lead1 =", cls(lead1), ",", cls(mon), index)
                square1 = {mymath.add_tuple(m, dif) for m in square}
                non_square1 = {mymath.add_tuple(m, dif) for m in non_square}
                if index == 2:
                    non_square1.remove(lcm)
                    dif1 = mymath.sub_tuple(lcm, lcms2[lcm])
                    rels_nsq = {mymath.add_tuple(m, dif1) for m in rels_lead_non_square[lcms2[lcm]]}
                    non_square1 ^= rels_nsq
                    new_square = {m for m in non_square1 if cls._is_square(m)}
                    square1 ^= new_square
                    non_square1 -= new_square
                    assert min(non_square1) > lcm
                non_square1 = cls.simplify_data(non_square1)
                new_square = {m1 for m1 in non_square1 if cls._is_square(m1)}
                square1 ^= new_square
                non_square1 -= new_square
                sq = cls._sqrt_zero(lead1, square1, non_square1, rels_lead_non_square, deg_max_rel)  # recursive
                if sq:
                    return sq

    @classmethod
    def _reduce(cls):
        """Reduce a square-nilpotent element. Return True if find one, otherwise False."""
        leads = [tuple(i + (i & 1) for i in m) for m in cls._rels if any(i > 1 for i in m)]
        leads.sort(key=cls.deg_mon)
        rels_lead_non_square = {}
        deg_max_rel = max(map(cls.deg_mon, cls._rels))
        for m in cls._rels:
            if not cls._is_square(m):
                dif = tuple(i & 1 for i in m)
                rel = {mymath.add_tuple(_m, dif) for _m in cls._rels[m]}
                rel.add(mymath.add_tuple(m, dif))
                mon = cls._first_non_square(rel)
                if mon:
                    rel.remove(mon)
                    rels_lead_non_square[mon] = rel
        for lead in leads:
            print("lead =", cls(lead))
            non_square = cls.simplify_data({lead})
            square = {_m for _m in non_square if cls._is_square(_m)}
            non_square -= square
            sq = cls._sqrt_zero(lead, square, non_square, rels_lead_non_square, deg_max_rel)
            if sq:
                root = {tuple(i // 2 for i in m1) for m1 in sq}
                root = cls.simplify_data(root)
                print("  new rel =", cls(root))
                cls.add_rel(root)
                return True
        return False

    @classmethod
    def reduce(cls):
        """Reduce all nilpotent elements."""
        while cls._reduce():
            pass
        cls.simplify_rels()

    # getters --------------------------
    @classmethod
    def is_admissible(cls, mon):
        """Determine if mon is in the basis."""
        return not any(mymath.le_tuple(m, mon) for m in cls._rels)

    @classmethod
    def get_rel_gens(cls):
        """Return a basis of relations."""
        rels, cls._rels = cls._rels, {}
        rel_gens = []
        try:
            for mon in sorted(rels, key=cls.deg_mon):
                rel_data = rels[mon] | {mon}
                rel_data = cls.simplify_data(rel_data)
                if rel_data:
                    rel_gens.append(rel_data)
                    cls.add_rel(rel_data)
            return rel_gens
        finally:
            cls._rels = rels

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = cls._gen_names.index(k)
        m = (0,) * i + (1,)
        return cls(m).simplify()

    @classmethod
    def _basis_mons_max(cls, deg_max, n_max):
        """Return an iterator of basis with length n_max + 1 and possibly trailing zeroes."""
        if n_max == -1:
            yield (), 0
            return
        for m, d in cls._basis_mons_max(deg_max, n_max - 1):
            for e in range((deg_max - d if d else deg_max) // cls._gen_degs[n_max] + 1):
                m1 = m + (e,)
                if any(map(mymath.le_tuple, cls._rels, itertools.repeat(m1))):
                    break
                else:
                    yield m1, d + cls._gen_degs[n_max] * e

    @classmethod
    def basis_mons_max(cls, deg_max):
        """Return an iterator of basis."""
        yield ()
        for n_max in range(len(cls._gen_degs)):
            for m, d in cls._basis_mons_max(deg_max, n_max - 1):
                for e in range(1, (deg_max - d if d else deg_max) // cls._gen_degs[n_max] + 1):
                    m1 = m + (e,)
                    if any(map(mymath.le_tuple, cls._rels, itertools.repeat(m1))):
                        break
                    else:
                        yield m1

    @classmethod
    def basis_max(cls, deg_max):
        return map(cls, cls.basis_mons_max(deg_max))

    @classmethod
    def present(cls, show_all=False):
        s = ", ".join(cls._gen_names)
        print(f"Generators: ${s}$.\\\\")
        print("Relations:\\\\")
        if show_all:
            for m in cls._rels:
                print(cls(m), "=", cls(cls._rels[m]))
        else:
            for data in cls.get_rel_gens():
                lead = min(data)
                print(f"${cls(lead)} = {cls(data - {lead})}$\\\\")

    @classmethod
    def ann(cls, x):
        """Return a basis for the ideal {y | xy=0}."""
        pass  # Todo: implement this


class SubRing:
    def __init__(self, d_max):
        self.d_max = d_max
        self.data = [[] for _ in range(d_max + 1)]

    def generate_non_comm(self, gens):
        if len(gens) == 0:
            raise ValueError("Need at least one generator")
        else:
            type_gen = type(gens[0])
            self.data[0].append([type_gen.unit(), type_gen.unit().get_mon()])
            for gen in gens:
                self.add_gen_non_comm(gen)

    def generate_comm(self, gens):
        if len(gens) == 0:
            raise ValueError("Need at least one generator")
        else:
            type_gen = type(gens[0])
            self.data[0].append([type_gen.unit(), type_gen.unit().get_mon()])
            for gen in gens:
                self.add_gen_comm(gen)

    def __str__(self):
        result = ""
        for d in range(self.d_max + 1):
            result += str(d) + ":\n"
            for base, _ in self.data[d]:
                result += str(base) + "\n"
        return result

    def add_gen_non_comm(self, gen):
        deg = gen.deg()
        if deg == 0:
            raise ValueError("the generator should have positive degree")
        for d in range(0, self.d_max - deg + 1):
            for base, _ in self.data[d]:
                new_base = base * gen
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

                new_base = gen * base  # non-commutative
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

    def add_gen_comm(self, gen):
        deg = gen.deg()
        if deg == 0:
            raise ValueError("the generator should have positive degree")
        for d in range(0, self.d_max - deg + 1):
            for base, _ in self.data[d]:
                new_base = base * gen
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

    def basis(self, deg):
        return (gen for gen, _ in self.data[deg])


# noinspection PyUnresolvedReferences,PyArgumentList
class QuoRing:
    """ This class must be inherited with another ring """
    d_max = None
    ideal = None  # type: linalg.GradedVectorSpaceMod2

    def __mul__(self, other):
        prod = super().__mul__(other)
        res = self.ideal.res(prod)
        return res

    @classmethod
    def init(cls, d_max):
        cls.d_max = d_max
        cls.ideal = linalg.GradedVectorSpaceMod2(d_max)

    @classmethod
    def add_relations(cls, rels, get_mon=None):
        for rel in rels:
            if rel:
                deg = rel.deg()
                for d in range(cls.d_max - deg + 1):
                    cls.ideal.add_vectors((super(QuoRing, a).__mul__(rel) for a in super().basis(d)), deg + d, get_mon)

    @classmethod
    def basis_mons(cls, deg):
        for m in super().basis_mons(deg):
            if m not in cls.ideal.get_mons(deg):
                yield m

    @classmethod
    def basis(cls, deg):

        return (cls(m) for m in cls.basis_mons(deg))


class FreeModule:
    type_ring = None

    def __init__(self, data: dict):
        if type(data) is dict:
            self.data = dict((mon, coeff) for mon, coeff in data.items() if coeff)
        else:
            raise TypeError

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        data = copy.deepcopy(self.data)
        for key in other.data:
            if key in data:
                data[key] += other.data[key]
            else:
                data[key] = other.data[key].copy()
        return type(self)(data)

    def __iadd__(self, other):
        data = self.data
        for key in other.data:
            if key in data:
                data[key] += other.data[key]
            else:
                data[key] = other.data[key].copy()
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __sub__(self, other):
        data = copy.deepcopy(self.data)
        for key in other.data:
            if key in data:
                data[key] -= other.data[key]
            else:
                data[key] = -other.data[key]
        return type(self)(data)

    def __isub__(self, other):
        data = self.data
        for key in other.data:
            if key in data:
                data[key] -= other.data[key]
            else:
                data[key] = -other.data[key]
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __rmul__(self, other):
        if type(other) is self.type_ring:
            data = {}
            for key in self.data:
                data[key] = other * self.data[key]
            return type(self)(data)
        else:
            return NotImplemented

    def __str__(self):
        result = ""
        for key in sorted(self.data):
            if result != "":
                result += "+"
            str_coeff = str(self.data[key])
            if '+' in str_coeff:
                result += "({}){}".format(str_coeff, key)
            elif str_coeff == "1":
                result += str(key)
            else:
                result += "{}{}".format(str_coeff, key)
        if result == "":
            result += "0"
        return result

    @classmethod
    def set_ring(cls, type_ring):
        cls.type_ring = type_ring

    @classmethod
    def gen(cls, key):
        """ return a generator of type FM """
        return cls({key: cls.type_ring.unit()})

    def is_zero(self) -> bool:
        return len(self.data) == 0

    def get_mon(self) -> tuple:
        key = max(self.data)
        return key, self.data[key].get_mon()

    def has_mon(self, mon: tuple) -> bool:
        if mon[0] in self.data:
            return mon[1] in self.data[mon[0]].data
        else:
            return False


class FreeModuleMod2:
    type_ring = None

    def __init__(self, data: Union[set, Tuple[tuple, Union[str, int]]]):
        if type(data) is set:
            self.data = data
        elif type(data) is tuple:  # monomial
            self.data = {data}  # type: Set[Tuple[tuple, Union[str, int]]]
        else:
            raise TypeError("{} can not initialize {}".format(data, type(self).__name__))

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        return type(self)(self.data ^ other.data)

    def __iadd__(self, other):
        self.data ^= other.data
        return self

    def __sub__(self, other):
        return type(self)(self.data ^ other.data)

    def __isub__(self, other):
        self.data ^= other.data
        return self

    def __rmul__(self, other):
        if type(other) is self.type_ring:
            data = set()
            for mon in self.data:
                prod = other * self.type_ring(mon[0])
                data.symmetric_difference_update((m, mon[1]) for m in prod.data)
            return type(self)(data)
        else:
            return NotImplemented

    def __str__(self):
        result = ""
        for m in self.data:
            if result != "":
                result += "+"
            str_coeff = str(self.type_ring(m[0]))
            if str_coeff == "1":
                result += str(m[1])
            else:
                result += "{}{}".format(str_coeff, m[1])
        if result == "":
            result += "0"
        return result

    @classmethod
    def set_ring(cls, type_ring):
        cls.type_ring = type_ring

    @classmethod
    def gen(cls, key):
        """ return a generator of type FM """
        return cls(((), key))

    @classmethod
    def zero(cls):
        return cls(set())

    def is_zero(self) -> bool:
        return len(self.data) == 0

    def get_mon(self) -> tuple:
        return max(self.data)

    def has_mon(self, mon: tuple) -> bool:
        return mon in self.data


def alg_may(n_max):
    R = AugAlgMod2.new_alg()
    B = {}
    K = {}
    # add generators
    for s in range(n_max):
        exec(f"K[{s}] = R.add_gen('h_{s}', {1})")
    for t in range(1, n_max + 1):
        for s in reversed(range(0, t)):
            if s == t - 1:
                B[(s, s + 1)] = K[s] * K[s]
            else:
                exec(f"B[({s}, {t})] = R.add_gen('b^{s}_{t-s}', 2)")
    for s in range(n_max - 2):
        exec(f"K[({s}, {s+1})] = R.add_gen('h_{s}(1)', 2)")
    for s in range(n_max - 4):
        exec(f"K[({s}, {s+1}, {s+2})] = R.add_gen('h_{s}(1,2)', 3)")
    for s in range(n_max - 4):
        exec(f"K[({s}, {s+1}, {s+3})] = R.add_gen('h_{s}(1,3)', 3)")

    # add relations
    for j in range(n_max + 1):
        for s in reversed(range(0, n_max - j + 1)):
            t = s + j
            rel = sum((B[(s, k)] * B[(k, t)] for k in range(s + 1, t)), R.zero())
            R.add_rel(rel)
            # print(s, t)
    for s in range(n_max):
        rel = K[s] * K[s] + B[(s, s + 1)]
        R.add_rel(rel)
        # print(s)
    for s in range(n_max - 2):
        rel = K[(s, s + 1)] * K[(s, s + 1)] + B[(s, s + 3)] * B[(s + 1, s + 2)] + \
              B[(s, s + 2)] * B[(s + 1, s + 3)]
        R.add_rel(rel)
        # print(s)
    R.reduce()
    R.present()
    return R, B, K


def test():
    pass


if __name__ == "__main__":
    pass

# 140, 248, 283, 415, 436, 612, 600
