from itertools import combinations
from algebras import BaseAlgebras as BA


# TODO commutativity at odd prime
class PolySingZ(BA.BasePolySingZ):
    # -- BasePolySingZ -------------
    @staticmethod
    def str_mon(mon: int) -> str:
        if mon >= 10 or mon < 0:
            return "x^{{{}}}".format(mon)
        elif mon > 1:
            return "x^{}".format(mon)
        elif mon == 1:
            return "x"
        else:
            return "1"

    def repr_(self, clsname="PolySingZ"):
        result = " + ".join(f"{c} * {self.repr_mon(m, clsname)}" for m, c in self._sorted_mons())
        return result if result else f"{clsname}.zero()"

    @staticmethod
    def repr_mon(mon: int, clsname="PolySingZ"):
        return f"{clsname}.gen({mon})"

    @staticmethod
    def unit_data():
        return {0: 1}

    # methods
    def __rmul__(self, other):
        return self * other

    def evaluation(self, x):
        """ compute f(x) """
        deg = self.deg()
        result = self.coeff(deg)
        for i in range(deg - 1, -1, -1):
            result = result * x + self.coeff(i)
        return result


class PolySingZ_trun(PolySingZ):
    deg_max = None

    # -- BasePolySingZ -------------
    @classmethod
    def gen(cls, exp: int = 1):
        if exp < 0:
            raise ValueError("negative exponent")
        return cls({exp: 1}) if exp <= cls.deg_max else cls.zero()

    def repr_(self, clsname="PolySingZ_trun"):
        return super().repr_(clsname)

    @staticmethod
    def repr_mon(mon: int, clsname="PolySingZ_trun"):
        return PolySingZ.repr_mon(mon, clsname)

    @classmethod
    def mul_data(cls, data1, data2):
        result = {}
        for m1, c1 in data1.items():
            for m2, c2 in data2.items():
                prod = m1 + m2
                if prod <= cls.deg_max:
                    if prod in result:
                        result[prod] += c1 * c2
                    else:
                        result[prod] = c1 * c2
        return result

    @classmethod
    def square_data(cls, data):
        result = {}
        for m, c in data.items():
            prod = m * 2
            if prod <= cls.deg_max:
                if prod in result:
                    result[prod] += c * c
                else:
                    result[prod] = c * c
        for item1, item2 in combinations(data.items(), 2):
            m1, c1 = item1
            m2, c2 = item2
            prod = m1 + m2
            if prod <= cls.deg_max:
                if prod in result:
                    result[prod] += c1 * c2 * 2
                else:
                    result[prod] = c1 * c2 * 2
        return result

    # methods
    def __truediv__(self, other: "PolySingZ_trun"):
        inv = other.inverse(self.deg_max)
        return self * inv

    @classmethod
    def set_deg_max(cls, d_max):
        cls.deg_max = d_max

    def composition(self, other):
        """ compute the composition of two polynomials """
        deg = self.deg()
        result = self.unit() * self.coeff(deg)
        for i in reversed(range(deg)):
            result = result * other + self.unit() * self.coeff(i)
        return result

    def inverse_composition(self):
        """Return the inverse in the sense of composition."""
        if 1 not in self.data or self.data[1] != 1 or 0 in self.data:
            raise ValueError("leading term is not x")
        data = {}
        y_powers = []
        prod = self.unit()
        for i in range(self.deg_max - 1):
            prod *= self
            y_powers.append(prod)
        data[0] = 0
        data[1] = 1
        for n in range(2, self.deg_max + 1):
            data[n] = -sum(data[i] * y_powers[i-1].coeff(n) for i in range(1, n))
        return type(self)(data)


class PolyAnyVarZ(BA.BasePolyAnyVar, BA.AlgebraZ):
    """
    This is for multi-variable polynomials over Z
    self.data is a dictionary of the form (monomial, coeff)
    where monomial is a tuple(dict) of the form (gen, exp)
    """
    dict_deg_gen = {}

    def _sorted_mons(self):
        """ for __str__ """
        gens = set(gen for mon in self.data for gen in dict(mon))
        return sorted(((tuple((gen, mon[gen]) if gen in mon else (gen, 0) for gen in sorted(gens)), coeff)
                       for mon, coeff in self._dict_items()), reverse=True)

    # methods -------------
    def _dict_items(self):
        """ turn mon into a dict """
        for mon, coeff in self.data.items():
            yield dict(mon), coeff


class PolyAnyVarModP(BA.BasePolyAnyVar, BA.AlgebraModP):
    """
    This is for multi-variable polynomials over Z/p
    self.data is a dictionary of the form (monomial, coeff)
    where monomial is a tuple(dict) of the form (gen, exp)
    """
    dict_deg_gen = {}

    # -- BasePolyModP -----------
    def frob(self):
        data = dict((tuple((gen, self.PRIME * exp) for gen, exp in mon), coeff) for mon, coeff in self.data.items())
        return type(self)(data)

    def _sorted_mons(self):
        """ for __str__ """
        gens = set(gen for mon in self.data for gen in dict(mon))
        return sorted(((tuple((gen, mon[gen]) if gen in mon else (gen, 0) for gen in sorted(gens)), coeff)
                       for mon, coeff in self._dict_items()), reverse=True)

    # methods -------------
    def _dict_items(self):
        """ turn mon into a dict """
        for mon, coeff in self.data.items():
            yield dict(mon), coeff


class PolyAnyVarMod2(BA.BasePolyAnyVar, BA.AlgebraMod2):
    """
    This is for multi-variable polynomials over Z/2
    self.data is a set of ((gen, exp),...) which is a tuplized dictionary
    """
    dict_deg_gen = {}

    # -- BasePolyMod2 -----------
    def square(self):
        data = set(tuple((gen, 2 * exp) for gen, exp in mon) for mon in self.data)
        return type(self)(data)

    def _sorted_mons(self):
        """ for __str__ """
        gens = set(gen for mon in self.data for gen in dict(mon))
        return sorted((tuple((gen, mon[gen]) if gen in mon else (gen, 0) for gen in sorted(gens))
                       for mon in self._dict_mons()), reverse=True)

    # methods -------------
    def _dict_mons(self):
        """ turn mon into a dict """
        for mon in self.data:
            yield dict(mon)


def poincare_series(gens_poly, gens_ext, d_max):
    """
    Given the degrees of generators,
    return the poincare series truncated by d_max
    """
    result_pro = PolySingZ.unit()
    for d in gens_poly:
        if d <= d_max:
            result_pro = result_pro.mul_trun((PolySingZ.unit() - PolySingZ.gen(d)).inverse(d_max), d_max)
    for d in gens_ext:
        if d <= d_max:
            result_pro = result_pro.mul_trun((PolySingZ.unit() + PolySingZ.gen(d)), d_max)
    return result_pro

# 372, 338, 318, 210, 183
