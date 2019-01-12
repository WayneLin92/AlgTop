"""
Warning: The degrees in this document are half of the topological degrees
"""
from .mymath import choose_mod2


# HBUZ constructors ----------------------------
def m_b(n):
    """ n >= 0 """
    return tuple(0 for _ in range(n)) + (1,)


def b(n):
    return HBUZ(tuple(0 for _ in range(n)) + (1,))


def b_rd(n):
    return b(n) * hz(-1)


def hz(n):
    """ Return [n] """
    return HBUZ((n,))


def psi_b(n):
    data = {(m_b(i), m_b(n - i)) for i in range(n + 1)}
    return HBUZT2(data)


# CoHBU constructors ----------------------------------
def m_c(n):
    """" c_0 = 1 """
    if n == 0:
        return ()
    return tuple(0 for _ in range(n-1)) + (1,)


def c(n):
    return CoHBU(m_c(n))


def psi_c(n):
    data = {(m_c(i), m_c(n - i)) for i in range(n + 1)}
    return CoHBUT2(data)


# Classes ---------------------------------------------
class CoHBU(algebras.PolyMod2):
    """ This is for polynomials of c_1, c_2, ... mod 2
        as a Hopf algebra """

    # virtual functions -----------------
    # noinspection PyShadowingNames
    @staticmethod
    def str_mon(m):
        str_result = ""
        index = 1
        for i in m:
            if i > 1 or i < 0:
                str_result += "c_{{{0}}}^{{{1}}}".format(index, i)
            elif i == 1:
                str_result += "c_{{{0}}}".format(index)
            index += 1
        if m == ():
            str_result = "1"
        return str_result

    @staticmethod
    def deg_mon(m):
        """ return the degree of m """
        s = 0
        d = 1
        for i in m:
            s += d * i
            d += 1
        return s

    @staticmethod
    def is_gen(m):
        """ return if m is a generator with its index """
        if m == m_c(len(m)):
            return True, len(m)
        else:
            return False, None

    @staticmethod
    def virschiebung(m):
        """
        Return the Virschiebung of the monomial.
        Return None if it is zero~
        """
        if len(m) % 2 == 1:
            return None
        else:
            result = [0] * (len(m) // 2)
        for i in range(len(m)):
            if i % 2 == 0:
                if m[i] != 0:
                    return None
            else:
                result[(i-1) // 2] = m[i]
        return tuple(result)

    def psi(self):
        result = CoHBUT2(set())
        for m in self.data:
            product = result.unit()
            index = 1
            for i in m:
                product = product * (psi_c(index) ** i)
                index += 1
            result += product
        return result

    def unit(self):
        return type(self)({(0,)})

    # dual ------------------
    @staticmethod
    def complexity(m):
        return sum(m)

    # noinspection PyUnusedLocal
    @staticmethod
    def pair(index_m, n):
        return 1 if len(n) <= 2 else 0

    @staticmethod
    def type_dual():
        return HBUZ


class CoHBUT2(algebras.PolyT2Mod2):
    """ Tensor product of two CoHBU """
    type_component = CoHBU


class HBUZ(algebras.PolyMod2):
    """
    This is for polynomials of b_0, b_1, b_2, ... mod 2
    as a Hopf algebra.
    The degree used here is one half of the topological degree
    including the degree of the Dyer-Lashof operations.
    """
    def inverse(self, d_max):
        list_homo = self.split_homo(d_max)
        if len(list_homo[0].data) != 1:
            raise ValueError("Not invertible")
        else:
            a0 = next(iter(list_homo[0].data))

        if len(a0) != 1:
            raise ValueError("Not invertible")
        else:
            list_result = [type(self)(set()) for _ in range(d_max + 1)]  # type, list[type(self)]
            list_result[0] = type(self)((-a0[0],))
            for d in range(1, d_max + 1):
                list_result[d] = sum((list_result[i] * list_homo[d - i] * type(self)((-a0[0],))
                                      for i in range(0, d)), type(self)(set()))
        return list_result

    def psi(self):
        """ coproduct """
        result = HBUZT2(set())
        for m in self.data:
            product = result.unit()
            index = 0
            for i in m:
                if index == 0:
                    product = product * HBUZT2({((i,), (i,))})
                else:
                    product = product * (psi_b(index) ** i)
                index += 1
            result += product
        return result

    def Q(self, s):
        """
        Dyer-Lashof operations.
        Warning: need s/2 <= NUM_CONJUGATES if b_0^{-1} is involved
        s here is half of the topological degree
        """
        global count_test
        count_test += 1

        sum_result = self.zero()
        for m in self.data:
            m_is_gen, index_m = self.is_gen(m)
            if m_is_gen:
                if s > index_m:
                    if index_m > 0:
                        for i in range(index_m+1):
                            if choose_mod2(s - index_m + i - 1, i) == 1:
                                sum_result += b(s+i) * b(index_m-i)
                    else:
                        if m[0] == 1:
                            sum_result += b(0) * b(s)
                        elif m[0] == -1:
                            sum_result += hz(-1) * chi_b[s]
                        elif m[0] == 0:
                            pass  # since s > 0
                elif s == index_m:
                    sum_result += HBUZ(m).square()
            else:
                mr = self.sqrt(m)
                if mr is not None:
                    if s % 2 == 0:
                        sum_result += HBUZ(mr).Q(s // 2).square()
                else:
                    m1, m2 = self.decompose(m)  # recursive
                    sum_result = sum((HBUZ(m1).Q(i) * HBUZ(m2).Q(s-i)
                                      for i in range(self.deg_mon(m1), s - self.deg_mon(m2) + 1)),
                                     sum_result)
        return sum_result

    def mQ(self, s):
        """
        Multiplicative Dyer-Lashof operations
        warning: need s <= NUM_CONJUGATES if b_0^{-1} is involved
        """
        global count_test
        count_test += 1

        sum_result = HBUZ(set())
        for m in self.data:
            m_is_gen, index_m = self.is_gen(m)
            if m_is_gen or index_m == 0:
                if index_m == 0:
                    if m[0] == 1 and s == 0:
                        sum_result += hz(1)
                    elif m[0] == 0 and s == 0:
                        sum_result += hz(0)
                    elif m[0] == -1:
                        sum_result += b(s)
                    elif m[0] == 2:
                        sum_result += hz(3) * b(s)
                    elif m[0] > 2:
                        if m[0] % 2 == 0:
                            mr = m[0] // 2
                            sum_result = sum((hz(mr).mQ(i).square() * hz(mr * mr).Q(s - 2 * i)
                                              for i in range(s // 2 + 1)), sum_result)
                        else:
                            sum_result = sum((hz(1) * hz(m[0]-1).mQ(i) * hz(m[0]-1).Q(s - i)
                                              for i in range(s + 1)), sum_result)
                    elif m[0] < -2:
                        sum_result = sum((b(i) % hz(-m[0]).mQ(s-i) for i in range(s+1)), sum_result)
            else:
                mr = self.sqrt(m)
                if mr is not None:
                    for mr1, mr2 in HBUZ(mr).psi().data:
                        sum_result = sum((HBUZ(mr1).mQ(i).square() * HBUZ(mr2).square_circ().Q(s-2*i)
                                          for i in range(self.deg_mon(mr1), s // 2 - self.deg_mon(mr2) + 1)),
                                         sum_result)
                else:
                    m1, m2 = self.decompose(m)
                    for m11, m12 in HBUZ(m1).psi().data:  # recursive
                        for m21, m22 in HBUZ(m2).psi().data:
                            sum_result = sum((sum((HBUZ(m11).mQ(i) * HBUZ(m21).mQ(j)
                                                   * (HBUZ(m12) % HBUZ(m22)).Q(s - i - j)
                                                   for i in range(self.deg_mon(m11),
                                                                  s - j - self.deg_mon(m12) - self.deg_mon(m22) + 1)),
                                                  HBUZ(set()))
                                              for j in range(self.deg_mon(m21), s + 1)), sum_result)
        return sum_result

    def __mod__(self, other):
        """
        circle product
        b_i\circ b_j Tested
        b_0^m\circ b_j Tested
        b_rd(j)\circ b_rd(j) Tested
        warning: need index_m, index_n <= NUM_CONJUGATES if b_0^{-1} is involved
        """
        global count_test
        count_test += 1

        sum_result = HBUZ(set())
        for m in self.data:
            m_is_gen, index_m = self.is_gen(m)
            for n in other.data:
                n_is_gen, index_n = self.is_gen(n)
                if m_is_gen and n_is_gen:
                    if index_m > 0 and index_n > 0:
                        if choose_mod2(index_m + index_n, index_n) == 1:
                            sum_result += b(index_m + index_n)
                    elif index_m == 0 and index_n == 0:
                        sum_result += hz(m[0] * n[0])
                    else:
                        if index_n == 0:  # switch m and n to make sure index_m = 0
                            mm, nn = n, m  # warning
                            index_n = index_m
                        else:
                            mm, nn = m, n
                        if mm[0] == 1:
                            sum_result += b(index_n)
                        elif mm[0] == -1:
                            sum_result += chi_b[index_n]
                        elif mm[0] == 0:
                            pass  # since index_n > 0 here
                elif index_m == 0 and index_n == 0:
                    sum_result += hz(m[0] * n[0])
                else:
                    if m_is_gen:  # switch m and n to make sure m is decomposable
                        mm, nn = n, m
                    else:
                        mm, nn = m, n
                    mr = self.sqrt(mm)
                    if mr is not None:
                        nrv = self.virschiebung(nn)
                        if nrv is not None:
                            sum_result += (HBUZ(mr) % HBUZ(nrv)).square()  # recursive
                    else:
                        m1, m2 = self.decompose(mm)
                        psi_n = HBUZ(nn).psi()
                        for n1, n2 in psi_n.data:
                            sum_result += (HBUZ(m1) % HBUZ(n1)) * (HBUZ(m2) % HBUZ(n2))  # recursive
        return sum_result

    def square_circ(self):
        return sum((HBUZ(m) % HBUZ(m) for m in self.data), HBUZ(set()))

    # virtual functions that overwrite
    @staticmethod
    def str_mon(m):
        str_result = ""
        index = 0
        for i in m:
            if i > 1 or i < 0:
                str_result += "b_{{{0}}}^{{{1}}}".format(index, i)
            elif i == 1:
                str_result += "b_{{{0}}}".format(index)
            index += 1
        if m == (0,):
            str_result = "1"
        return str_result

    @staticmethod
    def deg_mon(m):
        """ return the degree of m """
        s = 0
        d = 0
        for i in m:
            s += d * i
            d += 1
        return s

    @staticmethod
    def is_gen(m):
        """ If m is a generator return (True, index) otherwise (False, None) """
        if len(m) == 1:
            if m[0] in (1, 0, -1):
                return True, 0
            else:
                return False, None
        else:
            if m == m_b(len(m) - 1):
                return True, len(m) - 1
            else:
                return False, None

    @staticmethod
    def virschiebung(m):
        """
        Return the Virschiebung of the monomial.
        Return None if it is zero
        """
        if len(m) % 2 == 0:
            return None
        else:
            result = [0] * ((len(m) + 1) // 2)
        for i in range(len(m)):
            if i % 2 == 1:
                if m[i] != 0:
                    return None
            else:
                result[i // 2] = m[i]
        return tuple(result)

    def str_rd(self):
        """ print using b_rd """
        s_result = ""
        plus = False
        list_p = sorted(self.data, key=lambda n: (self.deg_mon(n), n), reverse=True)
        for m in list_p:
            if plus:
                s_result += "+"
            else:
                plus = True
            index = 0
            comp = sum(m)
            for i in m:
                if index == 0:
                    if comp > 1 or comp < 0:
                        s_result += "b_{{{0}}}^{{{1}}}".format(index, comp)
                    elif comp == 1:
                        s_result += "b_{{{0}}}".format(index)
                else:
                    if i > 1 or i < 0:
                        s_result += "\\olb_{{{0}}}^{{{1}}}".format(index, i)
                    elif i == 1:
                        s_result += "\\olb_{{{0}}}".format(index)
                index += 1
            if m == (0,):
                s_result += "1"
        if s_result == "":
            s_result += "0"
        return s_result

    def unit(self):
        return type(self)({(0,)})

    # dual ---------------------------------
    @staticmethod
    def complexity(m):
        return sum(m)

    # noinspection PyUnusedLocal
    @staticmethod
    def pair(index_m, n):
        return 1 if len(n) == 1 else 0

    @staticmethod
    def type_dual():
        return CoHBU


class HBUZT2(algebras.PolyT2Mod2):
    """ Tensor product of two HBUZ """
    # virtual functions
    def unit(self):
        return type(self)({((0,), (0,))})

    type_component = HBUZ


# Other functions ---------------------------------------------
def monomials_c(degree, index_max=None):
    """
    Iterator for looping over monomials of c_1, c_2, ..., c_{index_max} of the same degree
    in the form of tuples~
    """
    if index_max is None:
        index_max = degree
    if index_max == 1:
        yield degree,
    else:
        for i in range(degree // index_max, -1, -1):
            for t in monomials_c(degree - i * index_max, index_max - 1):
                if i > 0:
                    yield t + (0,) * (index_max - len(t) - 1) + (i,)
                else:
                    yield t


def print_rd(a):
    print(a.str_rd())


# Global variables ------------------------------
NUM_CONJUGATES = 20
chi_b = sum((b(n) for n in range(NUM_CONJUGATES + 1)), HBUZ(set())).inverse(NUM_CONJUGATES)


# Test ------------------------------------------
count_test = 0


def count():
    print("\ncount_bu =", count_test)


def hopf_ring_bu_test():
    if False:  # Test HBUZ.mQ(s)
        r = 3
        for s in range(10):
            print(s, hz(r + 1).mQ(s))
            print(s, sum((hz(1) * hz(r).mQ(i) * hz(r).Q(s - i) for i in range(s + 1)), HBUZ(set())))
        print(hz(6).mQ(11))
        # print(sum((hz(2).mQ(j) % hz(3).mQ(11-j) for j in range(11+1)), HBUZ(set())))
        print(b_rd(1).square().actQ(8).str_rd())
    if False:  # Test HBUZ.opQ(s)
        r = 5
        for s in range(10):
            print(s, hz(r+1).Q(s))
            print(s, sum((hz(r).Q(i)*hz(1).Q(s-i) for i in range(s+1)), HBUZ(set())))
        for s in range(5):
            print(s, b_rd(1).actQ(s).str_rd())
        print(b_rd(2).actQ(3).str_rd())
        print(b_rd(2).actQ(6).str_rd())
    if False:  # Test circle product
        for s in range(5):
            print((b_rd(1) % b_rd(s)).str_rd())
        print(hz(3) % hz(5))
        print(hz(2) % b(4))
        print(hz(3) % b(1))
    if False:
        # Performance
        N = 10
        print(chi_b[N].square_circ())
        print(count_test)
    if False:  # polynomial_test for item
        for m in monomials_c(10):
            print(CoHBU(m), algebras.pair(b(2) ** 5, CoHBU(m)))

    if True:
        n = 4
        x = b_rd(1) % b_rd(n)
        y = b_rd(1) % chi_b[n]
        print_rd(x)
        print_rd(y)
    algebras.count()
    count()


if __name__ == "__main__":
    hopf_ring_bu_test()
