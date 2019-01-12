""" Check if the formula for H_*BU<n> is correct """
from algebras.mymath import two_expansion
from algebras.BU import *
from algebras.operations import *
from algebras.polynomials import *


def is_stong_basis(mon):
    """ Determine a basis for A/(ASq^1+ASq^3) """
    if mon == ():
        return True
    m = [mon[i] - 2 * mon[i+1] for i in range(len(mon) - 1)] + [mon[-1]]

    index = len(m) - 1
    while index >= 0 and m[index] in (0, 2):
        index -= 1
    if index == -1:
        return True
    elif m[index] in (1, 3, 5):
        return False
    elif m[index] == 4:
        index -= 1
        while index >= 0 and m[index] == 0:
            index -= 1
        if index == -1:
            return True
        elif m[index] == 1:
            return False
        else:
            return True
    else:
        return True


def bbb(index_list):
    result = HBUZ({(1,)})
    for index in index_list:
        result = result % b_rd(index)
    return result


def mon_circ(mon):
    """ mon: dict -> list """
    n, m = list(mon.keys()), list(mon.values())
    k = [two_expansion(ms) for ms in m]
    r = [len(ks) for ks in k]
    m_all = sum(m)
    r_all = sum(r)
    n_new = [[n[s] - m_all + r_all + k[s][i] for i in range(len(k[s]))]
             for s in range(len(n))]
    n_new = sum(n_new, [])
    mon_new = {}
    for ns in n_new:
        if ns < 0:
            return None, None, None
        elif ns in mon_new.keys():
            mon_new[ns] += 1
        else:
            mon_new[ns] = 1
    conjugate = True if r_all % 2 == 1 else False

    needs_iterate = False
    for ms in mon_new.values():
        if ms > 1:
            needs_iterate = True
    if needs_iterate:
        mon_keys, conj, e = mon_circ(mon_new)  # recursive
        return mon_keys, conj != conjugate, e + m_all - r_all
    else:
        return list(mon_new.keys()), conjugate, m_all - r_all


poincare_d_max = 70


def test_mon_circ():
    mon = {2: 2, 4: 1}
    n, conj, e = mon_circ(mon)
    print(n, conj, e)
    print(sum(m * (2 ** n) for n, m in mon.items()))
    print((2 ** e) * sum(2 ** ns for ns in n))


def poincare_coh_bu(p):
    gen1 = [[d + 2 * p for m in Steenrod.basis_mons(d)
             if (m == () or m[-1] not in {1, 3}) and Steenrod.excess(m) < 2 * p
             and is_stong_basis(m)]
            for d in range(0, poincare_d_max + 1 - 2 * p)]
    gen1 = sum(gen1, [])  # A/ASq^1+ASq^3
    gen2 = [2 * i for i in range(poincare_d_max // 2 + 1) if bin(2 * i - 1).count("1") > p]
    gen = sorted(gen1 + gen2)
    print(gen1)
    print(gen2)
    print(gen)
    print(len(gen))
    f = poincare_series(gen, [], poincare_d_max)
    print(f)


def v(n):
    e = 0
    while n % 2 == 0:
        n = n // 2
        e += 1
    return e


def bar(g, d_max):
    """ Return list of 2^k(g+1) that are less than d_max """
    g = g + 1
    result = []
    while g <= d_max:
        result.append(g)
        g *= 2
    return result


def poincare_hom_bu(p):
    poly_gen = [2 * d for d in range(1, poincare_d_max // 2 + 1) if bin(d).count("1") >= p - v(d)]
    if p >= 3:
        gens = [[2 * d for d in range(1, poincare_d_max // 2 + 1) if bin(d).count("1") >= n - v(d)]
                for n in range(2, p)]
        degens = [[d for d in gens[n - 2] if bin(d).count("1") < n] for n in range(2, p)]
        ext_gen = []
        for n in range(2, p):
            ext_gen = sum((bar(e, poincare_d_max) for e in ext_gen), [])
            ext_gen = sum((bar(e, poincare_d_max) for e in ext_gen), [])
            ext_gen += sum((bar(e + 1, poincare_d_max) for e in degens[n - 2]), [])
    else:
        ext_gen = []

    print(poly_gen)
    print(ext_gen)

    f = poincare_series(poly_gen, ext_gen, poincare_d_max)
    print(f)


poincare_coh_bu(6)
poincare_hom_bu(6)


# test_mon_circ()
