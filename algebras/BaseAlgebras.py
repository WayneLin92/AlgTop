from abc import ABC, abstractmethod
from typing import Tuple, Set, Dict, Hashable, Any, Optional, Union
import operator
from itertools import product, combinations
from functools import reduce
from algebras import mymath
# todo: check the class methods
# todo: type is
# Todo: avoid creating new user-defined objects inside a class
# todo: try __slots__
# todo: consider unit degree


class Algebra(ABC):
    # methods -------------------
    def __repr__(self) -> str:
        return self.__str__()

    def __bool__(self) -> bool:
        return bool(self.data)

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        return type(self)(self.add_data(self.data, other.data))

    def __sub__(self, other):
        return type(self)(self.sub_data(self.data, other.data))

    def __mul__(self, other):
        return type(self)(self.mul_data(self.data, other.data))

    def __pow__(self, n: int):
        power = self.data
        pro = self.unit_data()
        while n:
            if n & 1:
                pro = self.mul_data(pro, power)
            n >>= 1
            if n:
                power = self.square_data(power)
        return type(self)(pro)

    def _repr_markdown_(self):
        return f"${self}$"

    def copy(self):
        return type(self)(self.data.copy())

    @classmethod
    def unit(cls):
        return cls(cls.unit_data())

    @classmethod
    def zero(cls):
        return cls(cls.zero_data())

    def square(self):
        return type(self)(self.square_data(self.data))

    def inverse(self, d_max) -> list:
        """Return 1/self, listed by degrees up to `d_max`."""
        list_homo = self.split_homo(d_max)
        if not list_homo or list_homo[0] != self.unit():
            raise ValueError("not monic")
        result = [self.unit()]
        for d in range(1, d_max + 1):
            term_d = -sum((result[i] * list_homo[d - i] for i in range(0, d)), self.zero())
            result.append(term_d)
        return result

    def deg(self) -> Optional[int]:
        """Return the deg of the polynomial. Return None if it is zero."""
        return max(map(self.deg_mon, self.data)) if self.data else None

    # abstract -----------------
    @abstractmethod
    def __init__(self, data):
        self.data = data  # type: Any
        raise NotImplementedError

    @abstractmethod
    def __str__(self) -> str: pass

    @abstractmethod
    def repr_(self, clsname: str) -> str:
        """Return the representation (functions as the actual `__repr__`)."""
        pass

    @staticmethod
    @abstractmethod
    def unit_data(): pass

    @staticmethod
    @abstractmethod
    def zero_data(): pass

    @staticmethod
    @abstractmethod
    def mul_mons(mon1, mon2):
        """Return the product of two monomials."""
        pass

    @staticmethod
    @abstractmethod
    def add_data(data1, data2):
        """Return the sum as data"""
        pass

    @staticmethod
    @abstractmethod
    def sub_data(data1, data2):
        """Return the sum as data"""
        pass

    @staticmethod
    @abstractmethod
    def mul_data(data1, data2):
        """Return product as data."""
        pass

    @staticmethod
    @abstractmethod
    def square_data(data):
        """Return the square of monomials."""
        pass

    @staticmethod
    @abstractmethod
    def str_mon(mon) -> str:
        """Return the str for the monomial."""
        pass

    @staticmethod
    @abstractmethod
    def repr_mon(clsname, mon) -> str:
        """Return the representation for the monomial."""
        pass

    @abstractmethod
    def _sorted_mons(self) -> list:
        """Sort the monomials for __str__()."""
        pass

    @abstractmethod
    def homo(self, d) -> "Algebra":
        """Return the degree d homogeneous part."""
        pass

    @abstractmethod
    def split_homo(self, d_max) -> list:
        """Return up to degree d homogeneous parts."""
        pass

    @staticmethod
    @abstractmethod
    def deg_mon(mon: Hashable):
        """Return the degree of mon."""
        pass


class AlgebraDict(Algebra, ABC):
    """ Algebra with coefficients. Use dict as data """
    # -- Algebra -----------
    def __str__(self):
        result = ""
        for m, c in self._sorted_mons():
            if result != "" and c > 0:
                result += "+"
            elif c < 0:
                result += "-"
            if c != 1 and c != -1:
                result += str(abs(c))
            result += self.str_mon(m)
        return result if result else "0"

    def __eq__(self, other):
        return {k: v for k, v in self.data.items() if v} == {k: v for k, v in other.data.items() if v}

    def __neg__(self):
        data = dict((m, -c) for m, c in self.data.items())
        return type(self)(data)

    def __iadd__(self, other):
        if type(other) is int:
            other = self.scalar(other)
        for mon, coeff in other.data.items():
            if mon in self.data:
                self.data[mon] += coeff
            else:
                self.data[mon] = coeff
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __isub__(self, other):
        if type(other) is int:
            other = self.scalar(other)
        for mon, coeff in other.data.items():
            if mon in self.data:
                self.data[mon] -= coeff
            else:
                self.data[mon] = -coeff
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __mul__(self, other):
        if type(other) is int:
            if other == 0:
                return self.zero()
            else:
                data = {}
                for mon in self.data:
                    data[mon] = (self.data[mon] * other)
                return type(self)(data)
        elif type(other) is type(self):
            return type(self)(self.mul_data(self.data, other.data))
        else:
            return NotImplemented

    def __imul__(self, other):
        if type(other) is int:
            if other == 0:
                self.data = {}
                return self
            else:
                for mon in self.data:
                    self.data[mon] *= other
                return self
        elif type(other) is type(self):
            return self * other
        else:
            return NotImplemented

    @classmethod
    def unit_data(cls):
        return {(): 1}

    @classmethod
    def zero_data(cls):
        return {}

    @classmethod
    def add_data(cls, data1, data2):
        result = data1.copy()
        for mon, coeff in data2.items():
            if mon in result:
                result[mon] += coeff
            else:
                result[mon] = coeff
        return result

    @classmethod
    def sub_data(cls, data1, data2):
        result = data1.copy()
        for mon, coeff in data2.items():
            if mon in result:
                result[mon] -= coeff
            else:
                result[mon] = -coeff
        return result

    @classmethod
    def mul_data(cls, data1, data2):
        result = {}
        for m1, c1 in data1.items():
            for m2, c2 in data2.items():
                prod = cls.mul_mons(m1, m2)
                if type(prod) is dict:
                    for m, c in prod:
                        if m in result:
                            result[m] += c * c1 * c2
                        else:
                            result[m] = c * c1 * c2
                else:
                    if prod in result:
                        result[prod] += c1 * c2
                    else:
                        result[prod] = c1 * c2
        return result

    @classmethod
    def square_data(cls, data):
        result = {}
        for m, c in data.items():
            prod = cls.mul_mons(m, m)
            if type(prod) is dict:
                for m_prod, c_prod in prod:
                    if m_prod in result:
                        result[m_prod] += c_prod * c * c
                    else:
                        result[m_prod] = c_prod * c * c
            else:
                if prod in result:
                    result[prod] += c * c
                else:
                    result[prod] = c * c
        for item1, item2 in combinations(data.items(), 2):
            m1, c1 = item1
            m2, c2 = item2
            prod = cls.mul_mons(m1, m2)
            if type(prod) is dict:
                for m_prod, c_prod in prod:
                    if m_prod in result:
                        result[m_prod] += c_prod * c1 * c2 * 2
                    else:
                        result[m_prod] = c_prod * c1 * c2 * 2
            else:
                if prod in result:
                    result[prod] += c1 * c2 * 2
                else:
                    result[prod] = c1 * c2 * 2
        return result

    # -- GradedAlgebra --------
    def homo(self, d):
        data = dict((m, c) for m, c in self.data.items() if self.deg_mon(m) == d)
        return type(self)(data)

    def split_homo(self, d_max):
        list_homo = [self.zero() for _ in range(d_max + 1)]
        for m, c in self.data.items():
            if self.deg_mon(m) <= d_max:
                if m in list_homo[self.deg_mon(m)].data:
                    list_homo[self.deg_mon(m)].data[m] += c
                else:
                    list_homo[self.deg_mon(m)].data[m] = c
        return list_homo

    # methods ----------------
    @classmethod
    def scalar(cls, n):  # assuming data is tuple
        return cls({(): n})

    def _sorted_mons(self) -> list:
        return sorted(self.data.items(), key=lambda item: (self.deg_mon(item[0]), item), reverse=True)


class BasePolyMulti(Algebra, ABC):
    """ class for multi-var polynomials """
    # -- Algebra --------------
    mul_mons = staticmethod(mymath.add_dict)

    @classmethod
    def str_mon(cls, mon: tuple):
        result = "".join(mymath.tex_pow(cls.str_gen(gen), exp) for gen, exp in mon)
        if result == "":
            result = "1"
        return result

    # -- GradedAlgebra
    @classmethod
    def deg_mon(cls, mon: tuple) -> int:
        return sum(cls.deg_gen(gen) * exp for gen, exp in mon)

    # abstract --------------
    @classmethod
    @abstractmethod
    def gen(cls, key: Hashable):
        """ return the key'th generator """
        pass

    @staticmethod
    @abstractmethod
    def deg_gen(key: Hashable) -> int:
        """ return the degree of the generator"""
        pass

    @staticmethod
    @abstractmethod
    def str_gen(key: Hashable) -> str:
        """ return the string for the generator"""
        pass


class BaseExterior(Algebra, ABC):
    """Base class for multi-variable exterior algebra."""
    # -- Algebra --------------
    @staticmethod
    def mul_mons(mon1: frozenset, mon2: frozenset) -> frozenset:
        return set() if mon1 & mon2 else mon1 | mon2

    @classmethod
    def str_mon(cls, mon: frozenset) -> str:
        result = "".join(map(cls.str_gen, sorted(mon)))
        return result if result else "1"

    # -- GradedAlgebra ----------------
    @classmethod
    def deg_mon(cls, mon: frozenset) -> int:
        return sum(map(cls.deg_gen, mon))

    # abstract --------------
    @classmethod
    @abstractmethod
    def gen(cls, key: Hashable):
        """ return the key'th generator """
        pass

    @staticmethod
    @abstractmethod
    def deg_gen(key: Hashable) -> int:
        """ return the degree of the generator"""
        pass

    @staticmethod
    @abstractmethod
    def str_gen(key: Hashable) -> str:
        """ return the string for the generator"""
        pass


class BasePolyAnyVar(BasePolyMulti, ABC):
    """ class for multi-var polynomials with any generator names """
    dict_deg_gen = None

    @classmethod
    def gen(cls, key, deg=1):
        if key in cls.dict_deg_gen and deg != cls.dict_deg_gen[key]:
            print("Warning: the degree of {} is changed from {} to {}".format(key, cls.dict_deg_gen[key], deg))
        cls.dict_deg_gen[key] = deg
        return cls(((key, 1),))

    @staticmethod
    def str_gen(key: str):
        return key

    @classmethod
    def deg_gen(cls, key):
        return cls.dict_deg_gen[key]


class Operations(ABC):
    def square(self):
        # noinspection PyUnresolvedReferences
        return self * self

    # abstract ------------
    @abstractmethod
    def simplify(self, degree: Optional[int] = None):
        """ Simplify the expression by adem relations and unstable relations """
        pass

    @staticmethod
    @abstractmethod
    def is_null(mon: tuple, degree=None) -> bool:
        """ determine if mon is zero for degree reason """
        pass

    @staticmethod
    @abstractmethod
    def is_admissible(mon: tuple) -> Tuple[bool, Optional[int]]:
        """
        Determine if mon is admissible and return also the
        place to apply the adem relations to
        """
        pass

    @staticmethod
    @abstractmethod
    def adem(mon: tuple, index: int) -> Set[tuple]:
        """ apply adem relation to mon[index], mon[index+1] and return a set of tuples """
        pass

    @classmethod
    @abstractmethod
    def gen(cls, n: int):
        pass


class HopfAlgebra(ABC):
    """ for Hopf algebras.
        Interface for pairing
    """
    # abstract ------------------
    @classmethod
    @abstractmethod
    def coprod_gen(cls, n: int): pass

    @abstractmethod
    def coprod(self): pass

    @staticmethod
    @abstractmethod
    def sqrt(mon: tuple) -> Optional[tuple]: pass

    @staticmethod
    @abstractmethod
    def virschiebung(mon: tuple): pass

    @staticmethod
    @abstractmethod
    def is_gen(mon: tuple) -> Tuple[bool, int]: pass

    @staticmethod
    @abstractmethod
    def decompose(mon: tuple) -> Tuple[tuple, tuple]: pass

    @staticmethod
    @abstractmethod
    def type_T2():
        pass


# Integers -----------------------------------
class AlgebraZ(AlgebraDict, ABC):
    """Class for algebra over Z."""
    def __init__(self, data: Union[dict, tuple]):
        if type(data) is tuple:
            self.data = {data: 1}
        elif type(data) is dict:
            self.data = dict(mon_coeff for mon_coeff in data.items() if mon_coeff[1])  # type: Dict[tuple, int]
        else:
            raise TypeError("{} can not initialize {}.".format(data, type(self).__name__))


class BasePolyZ(BasePolyMulti, AlgebraZ, ABC):
    pass


class BasePolySingZ(AlgebraDict, ABC):
    # -- AlgebraDict -------------
    def __init__(self, data: dict):
        if type(data) is dict:
            self.data = dict((mon, coeff) for mon, coeff in data.items() if coeff)  # type: Dict[int, int]
        else:
            raise TypeError("{} can not initialize {}.".format(data, type(self).__name__))

    def deg(self):
        if len(self.data) > 0:
            return max(self.data)
        else:
            return None

    @classmethod
    def unit(cls):
        return cls({0: 1})

    @classmethod
    def scalar(cls, n):
        return cls({0: n})

    def _sorted_mons(self):
        return sorted(self.data.items())

    def inverse(self, d_max):
        if 0 not in self.data or self.data[0] != 1:
            raise ValueError("not monic")
        data = {0: 1}
        for d in range(1, d_max + 1):
            data[d] = -sum(data[i] * self.coeff(d - i) for i in range(0, d))
        return type(self)(data)

    deg_mon = None

    # -- Algebra -----------------------
    @staticmethod
    def mul_mons(mon1: int, mon2: int) -> int:
        return mon1 + mon2

    # methods
    @classmethod
    def gen(cls, exp: int = 1):
        if exp < 0:
            raise ValueError("negative exponent")
        return cls({exp: 1})

    def coeff(self, n):
        return self.data[n] if n in self.data else 0


# Odd prime ---------------------------------
def _inv(a, p):
    for i in range(1, p):
        if a * i % p == 1:
            return i
    return None


class AlgebraModP(AlgebraDict, ABC):
    """ class for graded rings over F_p """
    PRIME = 2
    INV = [None, 1]

    def __init__(self, data: Union[dict, tuple]):
        if type(data) is tuple:
            self.data = {data: 1}
        elif type(data) is dict:
            self.data = dict((mon, coeff % self.PRIME)
                             for mon, coeff in data.items() if coeff % self.PRIME)  # type: Dict[tuple, int]
        else:
            raise TypeError("{} can not initialize {}.".format(data, type(self).__name__))

    def __iadd__(self, other):
        if type(other) is int:
            other = self.scalar(other)
        for mon, coeff in other.data.items():
            if mon in self.data:
                self.data[mon] += coeff
            else:
                self.data[mon] = coeff
        self.data = dict((mon, coeff % self.PRIME) for mon, coeff in self.data.items() if coeff % self.PRIME)
        return self

    def __isub__(self, other):
        if type(other) is int:
            other = self.scalar(other)
        for mon, coeff in other.data.items():
            if mon in self.data:
                self.data[mon] -= coeff
            else:
                self.data[mon] = -coeff
        self.data = dict((mon, coeff % self.PRIME) for mon, coeff in self.data.items() if coeff % self.PRIME)
        return self

    def __imul__(self, other):
        if type(other) is int:
            if other % self.PRIME == 0:
                self.data = {}
                return self
            else:
                for mon in self.data:
                    self.data[mon] *= other
                    self.data[mon] %= self.PRIME
                return self
        elif type(other) is type(self):
            return self * other
        else:
            raise NotImplementedError

    # methods ----------------
    @classmethod
    def set_prime(cls, p):
        cls.PRIME = p
        cls.INV = [_inv(i, p) for i in range(p)]

    @classmethod
    def get_prime(cls):
        return cls.PRIME

    def frob(self):
        return sum((type(self)((m, c)) ** self.PRIME for m, c in self.data.items()), self.zero())


class BasePolyModP(BasePolyMulti, AlgebraModP, ABC):
    pass


# Even prime ----------------------------------
class AlgebraMod2(Algebra, ABC):
    """ self.data is a set of monomials """
    def __init__(self, data: Union[set, tuple, frozenset]):
        if type(data) is set:
            self.data = data
        elif type(data) in (tuple, frozenset):  # monomial
            self.data = {data}  # type: Set[Union[tuple, frozenset]]
        else:
            raise TypeError(f"{data} of type {type(data)} can not initialize {type(self).__name__}.")

    # -- Algebra -----------
    def __str__(self):
        result = " + ".join(map(self.str_mon, self._sorted_mons()))
        return result if result else "0"

    def repr_(self, clsname):  # TODO: fix `type(self)` issue
        result = " + ".join(map(self.repr_mon, self._sorted_mons()))
        return result if result else f"{clsname}.zero()"

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (self.deg_mon(m), m), reverse=True)

    def __neg__(self):
        return self.copy()

    def __iadd__(self, other):
        self.data ^= other.data
        return self

    def __isub__(self, other):
        self.data ^= other.data
        return self

    @staticmethod
    def add_data(data1, data2):
        return data1 ^ data2

    @staticmethod
    def sub_data(data1, data2):
        return data1 ^ data2

    @classmethod
    def mul_data(cls, data1, data2):
        return reduce(operator.xor, (pro if type(pro := cls.mul_mons(m, n)) is set else {pro}
                                     for m, n in product(data1, data2)), set())

    @staticmethod
    def unit_data() -> set:
        return {()}

    @staticmethod
    def zero_data() -> set:
        return set()

    @classmethod
    def square_data(cls, data):
        """Warning: non-commutative algebra should overwrite this."""
        return reduce(operator.xor, (pro if type(pro := cls.mul_mons(m, m)) is set else {pro}
                                     for m in data), set())

    def homo(self, d):
        data = set(m for m in self.data if self.deg_mon(m) == d)
        return type(self)(data)

    def split_homo(self, d_max):
        list_homo = [self.zero() for _ in range(d_max + 1)]
        for m in self.data:
            if self.deg_mon(m) <= d_max:
                list_homo[self.deg_mon(m)].data.add(m)
        return list_homo

    def is_homo(self):
        prev_deg = None
        for m in self.data:
            if prev_deg is None:
                prev_deg = self.deg_mon(m)
            elif self.deg_mon(m) != prev_deg:
                return False
        return True


class AlgebraT2Mod2(AlgebraMod2, ABC):
    type_c0: AlgebraMod2 = None
    type_c1: AlgebraMod2 = None

    # -- AlgebraMod2 --------------
    def mul_mons(self, mon1: tuple, mon2: tuple):  # todo: return type
        prod0 = self.type_c0.mul_mons(mon1[0], mon2[0])
        prod1 = self.type_c1.mul_mons(mon1[1], mon2[1])
        if type(prod0) is tuple or type(prod0) is frozenset:  # assume both are monomials
            return prod0, prod1
        else:
            return {(m0, m1) for m0 in prod0 for m1 in prod1}

    def deg_mon(self, mon: tuple):
        """ return the degree of mon """
        deg0 = self.type_c0.deg_mon(mon[0])
        deg1 = self.type_c1.deg_mon(mon[1])
        return deg0 + deg1

    def str_mon(self, mon: tuple):
        str0 = self.type_c0.str_mon(mon[0])
        str1 = self.type_c1.str_mon(mon[1])
        return str0 + "\\otimes " + str1

    @classmethod
    def unit(cls):
        return cls({((), ())})

    # methods ----------------
    @classmethod
    def tensor(cls, a, b):
        assert type(a) is cls.type_c0 and type(b) is cls.type_c1
        return cls(set((m, n) for m in a.data for n in b.data))


class BasePolyMod2(BasePolyMulti, AlgebraMod2, ABC):
    pass


class BaseExteriorMod2(BaseExterior, AlgebraMod2, ABC):
    def __init__(self, data: Union[set, frozenset]):
        if type(data) is set:
            self.data = data
        elif type(data) is frozenset:  # monomial
            self.data = {data}  # type: Set[frozenset]
        else:
            raise TypeError("{} can not initialize {}.".format(data, type(self).__name__))

    def _sorted_mons(self) -> list:
        return sorted(self.data, key=lambda m: (self.deg_mon(m), tuple(m)), reverse=True)


class OperationsMod2(Operations, AlgebraMod2, ABC):
    # -- AlgebraMod2 ----------
    def __mul__(self, other):
        if not isinstance(other, OperationsMod2):
            return NotImplemented
        else:
            # noinspection PyUnresolvedReferences
            return super().__mul__(other).simplify()

    @staticmethod
    def mul_mons(mon1: tuple, mon2: tuple) -> tuple:
        return mon1 + mon2

    @staticmethod
    def deg_mon(mon: tuple):
        return sum(mon)

    # methods
    def simplify(self, degree: Optional[int] = None):
        s = self.data.copy()
        self.data = set()
        while len(s) > 0:
            m = s.pop()
            if self.is_null(m, degree):
                continue
            is_adm, index = self.is_admissible(m)
            if is_adm:
                self.data ^= {m}
            else:
                s ^= self.adem(m, index)
        return self


class HopfAlgWithDualMod2(HopfAlgebra, AlgebraMod2, ABC):
    # abstract --------------------
    @staticmethod
    @abstractmethod
    def complexity(mon: tuple) -> int: pass

    @staticmethod
    @abstractmethod
    def pair_gen(n: int, mon_dual: tuple): pass

    @staticmethod
    @abstractmethod
    def type_dual(): pass

    # methods ---------------------
    @classmethod
    def pair_mon(cls, mon1: tuple, mon2: tuple):
        type1, type2 = cls, cls.type_dual()
        if type1.deg_mon(mon1) != type2.deg_mon(mon2):
            return 0
        if type2.complexity(mon2) < type1.complexity(mon1):
            mon1, mon2 = mon2, mon1
            type1, type2 = type2, type1
        result = 0
        mon1_is_gen, gen1 = type1.is_gen(mon1)
        if mon1_is_gen:
            result = type1.pair_gen(gen1, mon2)
        else:
            root1 = type1.sqrt(mon1)
            if root1 is not None:
                virch2 = type2.virschiebung(mon2)
                if virch2 is not None:  # recursive
                    result += type1.pair_mon(root1, virch2)
            else:
                m11, m12 = type1.decompose(mon1)
                psi2 = type2(mon2).coprod()
                for m21, m22 in psi2.data:  # recursive
                    result += type1.pair_mon(m11, m21) * type1.pair_mon(m12, m22)
        return result % 2

    def pair(self, other):
        if self.type_dual() is not type(other):
            raise TypeError("{} and {} can not be paired.".format(self, other))
        result = sum(self.pair_mon(mon1, mon2) for mon1 in self.data for mon2 in other.data)
        return result % 2


# Monitor ---------------------------------------
class Monitor:
    dict_num_function_calls = {}

    @classmethod
    def count(cls, key):
        if key not in cls.dict_num_function_calls:
            cls.dict_num_function_calls[key] = 1
        else:
            cls.dict_num_function_calls[key] += 1

    @classmethod
    def present(cls):
        frequent_functions = sorted(cls.dict_num_function_calls.items(), key=lambda item: -item[1])
        print("\nnumber of function calls = {}".format(frequent_functions[:30]))


# 670, 615, 627, 643, 660, 693, 756, 764, 845, 851, 888, 884, 868, 885, 874, 852
