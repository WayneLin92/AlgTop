"""Classes that can be constructed by other classes."""
# todo: create class AugModuleMod2
# Todo: construct AugAlgMod2 from other algebras
import copy
import itertools
import operator
import heapq
from typing import Union, Set, Tuple, List, Dict
from algebras import BaseAlgebras as BA, linalg, mymath, myerror


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
    def new_alg(unit_deg=None) -> "AugAlgMod2":
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
    def add_gens(cls, iterable_names):
        """Add generators."""
        cls._gen_names += iterable_names

    @classmethod
    def add_rel(cls, rel: "AugAlgMod2"):
        """Add a relation."""
        hq = [(rel.deg(), rel.data)]
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
                            new_rel = set(map(mymath.add_tuple, r, itertools.repeat(dif)))
                            heapq.heappush(hq, (cls.deg_mon(lcm), new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls._rels[m_redundant]
                cls._rels[m] = r - {m}
                print(cls(m))

    @classmethod
    def simplify_rels(cls):
        for m in cls._rels:
            cls._rels[m] = cls.simplify_data(cls._rels[m])

    @classmethod
    def simplify_data(cls, data: set):
        """Simplify the data by relations."""
        s = data.copy()
        result = set()
        while len(s) > 0:
            mon = s.pop()
            for m in cls._rels:
                if mymath.le_tuple(m, mon):
                    q, r = mymath.div_mod_tuple(mon, m)
                    mon1 = set(mymath.add_tuple(tuple(map(operator.mul, m1, itertools.repeat(q))), r)
                               for m1 in cls._rels[m])
                    s ^= mon1
                    break
            else:
                result ^= {mon}
        return result

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    # getters --------------------------
    @classmethod
    def get_rel_gens(cls):
        rels, cls._rels = cls._rels, {}
        rel_gens = []
        try:
            for mon in sorted(rels, key=cls.deg_mon):
                rel_data = rels[mon] | {mon}
                rel_data = cls.simplify_data(rel_data)
                if rel_data:
                    rel_gens.append(rel_data)
                    cls._add_rel(rel_data)
            return rel_gens
        finally:
            cls._rels = rels

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = 0
        while i < len(cls._gen_names) and cls._gen_names[i] != k:
            i += 1
        if i == len(cls._gen_names):
            raise myerror.MyKeyError("Generator not found.")
        m = (0,) * (i - 1) + (1,)
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
    def present(cls):
        s = ", ".join(cls._gen_names)
        print(f"Generators: ${s}$.\\\\")
        print("Relations:\\\\")
        for data in cls.get_rel_gens():
            print(f"${cls(data)}=0$\\\\")

    @classmethod
    def nil(cls):
        """Return a basis for the nilradical."""
        leading_mons = sorted(cls._rels)
        leading_masks = [tuple(map(int, map(bool, m))) for m in leading_mons]
        for m in leading_mons:
            mask_odd = map(operator.mod, m, itertools.repeat(2))
            m1 = tuple(map(operator.add, m, mask_odd))
            square = {m1}

    @classmethod
    def ann(cls, x):
        """Return a basis for the ideal {y | xy=0}."""
        pass  # TODO: implement this


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


def test():
    R = AugAlgMod2.new_alg()
    B = {}
    n_max = 7
    for t in range(1, n_max + 1):
        for s in reversed(range(0, t)):
            exec(f"B[({s}, {t})] = R.add_gen('B_{{{s}{t}}}', {2 ** (t + 1) - 2 ** (s + 1)})")
    for j in range(2, n_max + 1):
        for s in reversed(range(0, n_max - j + 1)):
            t = s + j
            rel = sum((B[(s, k)] * B[(k, t)] for k in range(s + 1, t)), R.zero())
            print(s, t)
            R.add_rel(rel)
        R.simplify_rels()
    BA.Monitor.present()
    return R, B


if __name__ == "__main__":
    import timeit
    print(timeit.timeit("test()", "from __main__ import test", number=100))

# 140, 248, 283, 415, 436
