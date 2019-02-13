from algebras import BaseAlgebras as BA, myerror


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

    # methods
    def mul_trun(self, other, d_max):
        data = {}
        for m1, c1 in self.data.items():
            for m2, c2 in other.data.items():
                prod = m1 + m2
                if prod <= d_max:
                    if prod in data:
                        data[prod] += c1 * c2
                    else:
                        data[prod] = c1 * c2
        return type(self)(data)

    def inv_function(self, d_max):
        """ inverse in the sense of composition """
        if 1 not in self.data or self.data[1] != 1 or 0 in self.data:
            raise myerror.MyValueError("require the leading term to be x")
        data = {}
        y_powers = []
        prod = self.unit()
        for i in range(d_max-1):
            prod = prod.mul_trun(self, d_max)
            y_powers.append(prod)
        data[0] = 0
        data[1] = 1
        for n in range(2, d_max + 1):
            data[n] = -sum(data[i] * y_powers[i-1].coeff(n) for i in range(1, n))
        return type(self)(data)

    def composition(self, other, d_max):
        """ compute the composition of two polynomials """
        deg = self.deg()
        result = self.unit() * self.coeff(deg)
        for i in range(deg - 1, -1, -1):
            result = result.mul_trun(other, d_max) + self.coeff(i)
        return result

    def eval(self, x):
        """ compute f(x) """
        deg = self.deg()
        result = self.coeff(deg)
        for i in range(deg - 1, -1, -1):
            result = result * x + self.coeff(i)
        return result

    def pow(self, n, d_max):
        power = self
        pro = self.unit()
        while n > 0:
            if n % 2 == 1:
                pro = pro.mul_trun(power, d_max)
            n //= 2
            power = power.square_trun(d_max)
        return pro

    def square_trun(self, d_max):
        data = {}
        for m, c in self.data.items():
            prod = self.mul_mons(m, m)
            if prod <= d_max:
                if prod in data:
                    data[prod] += c * c
                else:
                    data[prod] = c * c
        list_data = list(self.data.items())
        for i in range(len(list_data)):
            for j in range(i + 1, len(list_data)):
                m1, c1 = list_data[i]
                m2, c2 = list_data[j]
                prod = self.mul_mons(m1, m2)
                if prod <= d_max:
                    if prod in data:
                        data[prod] += c1 * c2 * 2
                    else:
                        data[prod] = c1 * c2 * 2
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
