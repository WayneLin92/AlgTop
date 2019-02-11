""" classes that are constructed by other classes """
import copy
import itertools
import operator
from typing import Union, Set, Tuple, List, Dict
from algebras import BaseAlgebras as BA, linalg, mymath


class AugAlgMod2(BA.AlgebraMod2):
    """A factory for augmented algebras over F_2.

    AugAlgMod2.new_alg() creates a new augmented algebra with
    its own generators and relations."""

    _gen_names = None  # type: List[str]
    _gen_degs = None  # type: List[int]
    _null_mons = None  # type: List[tuple]
    _rels = None  # type: Dict[tuple, set]
    _name_index = 0

    # ----- AlgebraMod2 -------------
    @classmethod
    def str_mon(cls, mon: tuple):
        if mon:
            return "".join(map(lambda ie: f"{cls._gen_names[ie[0]]}{mymath.tex_exponent(ie[1])}", enumerate(mon)))
        else:
            return "1"

    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple):
        if len(mon1) < len(mon2):
            return tuple(itertools.chain(map(operator.add, mon1, mon2), mon2[len(mon1):]))
        else:
            return tuple(itertools.chain(map(operator.add, mon1, mon2), mon1[len(mon2):]))

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(lambda ie: cls._gen_degs[ie[0]] * ie[1], enumerate(mon)))

    # methods --------------------
    @staticmethod
    def new_alg():
        """Return a dynamically created subclass of AugAlgMod2."""
        cls = AugAlgMod2
        class_name = f"AugAlgMod2_{cls._name_index}"
        cls._name_index += 1
        dct = {'_gen_names': [], '_gen_degs': [], '_null_mons': [], '_rels': {}}
        return type(class_name, (cls,), dct)

    @classmethod
    def add_gen(cls, k: str, deg):
        """Add a new generator and return it."""
        cls._gen_names.append(k)
        cls._gen_degs.append(deg)
        m = (0,) * (len(cls._gen_names) - 1) + (1,)
        return cls(m)

    @classmethod
    def add_gens(cls, iterable_names):
        """Add generators."""
        cls._gen_names += iterable_names

    @classmethod
    def add_rel(cls, rel: "AugAlgMod2"):
        """Add a relation. Assert rel is simplified."""
        if len(rel.data) == 1:
            for m in rel.data:
                cls._null_mons.append(m)

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = 0
        while i < len(cls._gen_names) and cls._gen_names[i] != k:
            i += 1
        if i == len(cls._gen_names):
            raise BA.MyKeyError("Generator not found.")
        m = (0,) * (i - 1) + (1,)
        return cls(m)

    @classmethod
    def simplify_data(cls, data: set):
        s = data.copy()
        result = set()
        while len(s) > 0:
            mon = s.pop()
            for m in cls._null_mons:
                if mymath.leq_tuple(m, mon):
                    continue
            for m in cls._rels:
                if mymath.leq_tuple(m, mon):
                    q, r = mymath.div_mod_tuple(mon, m)

            is_adm, index = self.is_admissible(mon)
            if is_adm:
                self.data ^= {mon}
            else:
                s ^= self.adem(mon, index)
        return self


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


# 140, 248, 283
