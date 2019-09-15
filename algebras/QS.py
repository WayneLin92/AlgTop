"""
To be maintained.
The homology of $QS^0$.
"""
from .mymath import choose_mod2
from . import operations as op
from . import algebras


# constructors ----------------------------------
def hz(n):
    """ Return [n] """
    return HQS.hz(n)


# Classes ---------------------------------------------
class HQS(op.DyerLashofX):
    """
    This is for the homology of QS^0
    Data is a set of tuples
    Each tuple is actually a dict
    """
    def __rmul__(self, other):
        if type(other) is HQS:
            return super(HQS, self).__mul__(other)
        elif type(other) is op.DyerLashof:
            return sum((self.actQ(m) for m in other.data), self.zero())
        elif type(other) is op.Steenrod:
            return sum((self.actSq(m) for m in other.data), self.zero())
        elif type(other) is op.AR:
            return sum((self.actSq(n).actQ(m) for m, n in other.data), self.zero())
        else:
            return NotImplemented

    @classmethod
    def set_x_deg(cls, x_deg):
        raise NotImplementedError

    @classmethod
    def hz(cls, n):
        """ return [n] """
        return cls((((), n),))

    @classmethod
    def str_mon(cls, mon):
        result = ""
        for k, v in mon:
            if v > 1:
                if k != ():
                    if v >= 10:
                        result += "({}[1])^{{{}}}".format(op.DyerLashof(k), v)
                    else:
                        result += "({}[1])^{}".format(op.DyerLashof(k), v)
                else:
                    result += "[{}]".format(v)
            elif v == 1:
                if k != ():
                    result += "{}[1]".format(op.DyerLashof(k))
                else:
                    result += "[1]"
        if result == "":
            result = "[0]"
        return result

    def actSq(self, r):
        count()

        if type(r) is tuple:
            result = self
            for i in range(len(r)):  # Note that Sq^I is the dual action
                result = result.actSq(r[i])
            return result

        result = self.zero()
        for m in self.data:
            if len(m) == 1 and m[0][1] == 1:  # #########
                q = (op.Sq(r) * op.DyerLashof(m[0][0])).truncate(0)
                result += HQS(set(((m, 1),) for m, n in q.data))
            elif len(m) == 0:
                if r == 0:
                    result += self.unit()
            else:
                mr = self.sqrt(m)
                if mr is not None:
                    if r % 2 == 0:
                        result += HQS(mr).actSq(r // 2).square()
                else:
                    m1, m2 = self.decompose(m)
                    result = sum((HQS(m1).actSq(i) * HQS(m2).actSq(r - i)
                                  for i in range(self.deg_mon(m1) // 2, r - self.deg_mon(m2) // 2 + 1)), result)
        return result.simplify()

    def __mod__(self, other):
        """
        circle product
        b_i\circ b_j Tested
        b_0^m\circ b_j Tested
        b_rd(j)\circ b_rd(j) Tested
        warning: need index_m, index_n <= NUM_CONJUGATES if b_0^{-1} is involved
        """
        count()

        sum_result = self.zero()
        for m in self.data:
            m_is_gen, index_m = self.is_gen(m)
            for n in other.data:
                n_is_gen, index_n = self.is_gen(n)
                if m_is_gen and n_is_gen:
                    if len(index_n) < len(index_m):
                        index_m, index_n = index_n, index_m
                    if index_m == ():
                        sum_result += HQS(((index_n, 1),))
                    else:
                        x, y = HQS(((index_m[1:], 1),)), HQS(((index_n, 1),))
                        sum_result = sum((op.Q(index_m[0] + i) * (x % y.actSq(i))
                                          for i in range(self.deg_mon(((index_n, 1),)) // 2 + 1)), sum_result)
                else:
                    if m_is_gen:  # switch m and n to make sure m is decomposable
                        m, n = n, m
                    mr = self.sqrt(m)
                    if mr is not None:
                        nrv = self.virschiebung(n)
                        if nrv is not None:
                            sum_result += (HQS(mr) % HQS(nrv)).square()  # recursive
                    else:
                        m1, m2 = self.decompose(m)
                        psi_n = HQS(n).psi()
                        for n1, n2 in psi_n.data:
                            sum_result += (HQS(m1) % HQS(n1)) * (HQS(m2) % HQS(n2))  # recursive
        return sum_result

    @staticmethod
    def decompose(m):  # ######
        """
        Decompose a monomial into two monomials
        m1 is a generator or a power of a generator
        m must be decomposable
        """
        return (m[0:1], m[1:]) if len(m) > 1 else (((m[0][0], 1),), ((m[0][0], m[0][1] - 1),))

    @staticmethod
    def is_gen(m):
        """ return if m is a generator with its index """
        if len(m) == 1 and m[0][1] == 1:
            return True, m[0][0]
        else:
            return None, None

    @staticmethod
    def virschiebung(m):
        """
        Return the Virschiebung of the monomial.
        Return None if it is zero
        """
        for k, v in m:
            for i in k:
                if i % 2 == 1:
                    return None
        return tuple((tuple(i // 2 for i in k), v) for k, v in m)

    def psi(self):
        """ coproduct """
        count()

        result = HQST2(set())
        for m in self.data:
            product = result.unit()
            for k, v in m:
                if len(k) == 0:
                    product = product * HQST2((((k, v),), ((k, v),)))
                else:
                    coprod = HQS(((k[1:], v),)).psi()  # recursive
                    t = sum((HQST2.tensor(HQS(x1).actQ(i), HQS(x2).actQ(k[0] - i))
                            for x1, x2 in coprod.data for i in range(HQS.deg_mon(x1), k[0] - HQS.deg_mon(x2) + 1)),
                            result.zero())
                    product = product * (t ** v)
            result += product
        return result


class HQST2(algebras.PolyT2Mod2):
    """ Tensor product of two CoHBU """
    type_component = HQS


## Test
count_test = 0


def count():
    global count_test
    count_test += 1


def print_count():
    print("\ncount_qs =", count_test)


def qs_test():
    x = HQS.hz(3)
    y = op.Q(20) * (op.Q(8) * x)
    print(y)
    print(op.Q(0) * x)


def qs_virschiebung_test():
    x = HQS.hz(1)
    y = op.Q(20) * (op.Q(8) * x) + op.Q(21) * op.Q(8) * x
    print(y)
    for m in y.data:
        print(HQS.virschiebung(m))
    print("supposed to be (((9, 5), 1),) None None\n")


def qs_psi_test():
    x = HQS.hz(1)
    y = op.Q(6) * op.Q(3) * op.Q(2) * x
    print(y)
    print(y.psi())


def qs_mod_test():
    i1, j1 = 10, 9
    x = op.Q(i1) * HQS.hz(1)
    y = op.Q(j1) * HQS.hz(1)
    hz1 = HQS.hz(1)
    print(hz1 % hz1)
    print(hz1 % y)
    print('\n')
    print(x % y)

    k_max = 10
    r1 = op.sQ(j1, k_max) * op.sQ(i1, k_max)
    r2 = op.sQ(i1, k_max) * op.sQ(j1, k_max)
    print(r1 == r2)
    print(r1 * op.Q(4) * hz1)
    print(r2 * op.Q(4) * hz1)
    print(r1 * op.Q(4))
    print(r2 * op.Q(4))


if __name__ == "__main__":
    # qs_test()
    # bu_virschiebung_test()
    # qs_psi_test()
    qs_mod_test()
    print_count()
