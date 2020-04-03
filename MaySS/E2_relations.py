from itertools import permutations, combinations, product
from algebras.mymath import prod_algs
from algebras.groebner import GbAlgMod2


def key(m):
    return [-i for i in m]


def def_gens(R: GbAlgMod2):
    """Bind some local variables to the generators of `R`."""
    for name in R.gen_names:
        translation = str.maketrans("(", "_", " _){},")
        var_name = name.translate(translation)
        var_dict = globals()
        var_dict[var_name] = R.gen(name)


def decompose(S, T):
    """Decompose h_{S,T} into indecomposables."""
    left, right = min(S), max(T)
    i = left
    for j in range(left + 2, right + 2, 2):
        S1, T1 = S & set(range(i, j)), T & set(range(i, j))
        if len(S1) == len(T1):
            yield frozenset(S1), frozenset(T1)
            i = j


def h(S, T):
    R"""Return $h_i(\tilde S)$"""
    result = A.unit()
    if not S:
        return result
    left, right = min(S), max(T)
    if left > min(T):
        return A.zero()
    i = left
    for j in range(left + 2, right + 2):
        S1, T1 = S & set(range(i, j)), T & set(range(i, j))
        if len(S1) == len(T1):
            S1 = sorted(S1)
            iS1 = S1[0]
            tS = [S1[j] - iS1 for j in range(1, len(S1))]
            result *= A.gen(f"h_{iS1}({', '.join(map(str, tS))})") if tS else A.gen(f"h_{iS1}")
            i = j
        elif len(S1) < len(T1):
            return A.zero()
    return result


def b(S, T):
    """Return $b_{S, T}$"""
    if type(S) is int:
        if T == S + 1:
            return A.gen(f"h_{S}") ** 2
        else:
            return A.gen(f"b_{{{S}{T}}}")
    S, T = sorted(S), sorted(T)
    result = A.zero()
    for T1 in permutations(T):
        if all(t - s > 0 for s, t in zip(S, T1)):
            prod = A.unit()
            for s, t in zip(S, T1):
                prod *= b(s, t)
            result += prod
    return result


def h_interior(S, T):
    return h(S - {min(S)}, T - {max(T)})


def h_wrap(S, T):
    return h(S | {min(S) - 1}, T | {max(T) + 1})


def relation1(i: int, j: int):
    return sum((b(i, k) * b(k, j)for k in range(i + 1, j)), A.zero())


def relation2(S1: set, T1: set, S2: set, T2: set):
    return h(S1, T1) * h(S2, T2)


def relation3A(j: int, S: set, T: set):
    return sum((b(s, j) * h(S - {s}, T | {s}) for s in S if s < j), A.zero())


def relation3B(i: int, S: set, T: set):
    return sum((b(i, t) * h(S | {t}, T - {t}) for t in T if t > i), A.zero())


def relation4A(S1: set, T1: set, S2: set, T2: set):
    S1p = S1 - (S2 | T2)
    T1p = T1 - (S2 | T2)
    S1pp, T1pp = S1 - S1p, T1 - T1p
    return sum((h(S1pp | setI, T1pp - setI) * h(S1p | S2 - setI, T1p | T2 | setI)
                for setI in map(set, combinations(T1pp & S2, (len(T1pp) - len(S1pp)) // 2))),
               h(S1, T1) * h(S2, T2))


def relation4B(S1: set, T1: set, S2: set, T2: set):
    S2p = S2 - (S1 | T1)
    T2p = T2 - (S1 | T1)
    S2pp, T2pp = S2 - S2p, T2 - T2p
    return sum((h(S2pp - setI, T2pp | setI) * h(S2p | S1 | setI, T2p | T1 - setI)
                for setI in map(set, combinations(T1 & S2pp, (len(S2pp) - len(T2pp)) // 2))),
               h(S1, T1) * h(S2, T2))


def relation5(S1: set, T1: set, S2: set, T2: set):
    S1p = S1 - (S2 | T2)
    T1p = T1 - (S2 | T2)
    S1pp, T1pp = S1 - S1p, T1 - T1p
    S2p = S2 - (S1 | T1)
    T2p = T2 - (S1 | T1)
    S2pp, T2pp = S2 - S2p, T2 - T2p
    return sum((h(S1p - setI, T1p | setI) * b(S1pp | setI, T2pp | setJ) * h(S2p | setJ, T2p - setJ)
                for setI in map(set, combinations(S1p, (len(S1p) - len(T1p)) // 2))
                for setJ in map(set, combinations(T2p, (len(T2p) - len(S2p)) // 2))),
               h(S1, T1) * h(S2, T2))


def relation6(S1: set, T1: set, S2: set, T2: set):
    S1p, T1p = S1 - {min(S1)}, T1 - {max(T1)}
    S2p, T2p = S2 - {min(S2)}, T2 - {max(T2)}
    indecomposables1 = set(decompose(S1p, T1p))
    indecomposables2 = set(decompose(S2p, T2p))
    x = prod_algs((h(S, T) for S, T in (indecomposables2 - indecomposables1)), A.unit())
    y = prod_algs((h(S, T) for S, T in (indecomposables1 - indecomposables2)), A.unit())
    return x * h(S1, T1) + y * h(S2, T2)


def relation6AA(S1: set, T1: set, S2: set, T2: set, S3: set, T3: set):
    S1p = S1 - (S2 | T2)
    T1p = T1 - (S2 | T2)
    S1pp, T1pp = S1 - S1p, T1 - T1p
    yield (S1, T1), (S2 | S3, T2 | T3)
    for setI in map(set, combinations(T1pp & S2, (len(T1pp) - len(S1pp)) // 2)):
        yield (S1pp | setI | S3, T1pp - setI | T3), (S1p | S2 - setI, T1p | T2 | setI)


def relation6AB(S1: set, T1: set, S2: set, T2: set, S3: set, T3: set):
    S1p = S1 - (S2 | T2)
    T1p = T1 - (S2 | T2)
    S1pp, T1pp = S1 - S1p, T1 - T1p
    yield (S2, T2), (S1 | S3, T1 | T3)
    for setI in map(set, combinations(T1pp & S2, (len(T1pp) - len(S1pp)) // 2)):
        yield (S1pp | setI | S3, T1pp - setI | T3), (S1p | S2 - setI, T1p | T2 | setI)


def relation6BA(S1: set, T1: set, S2: set, T2: set, S3: set, T3: set):
    S2p = S2 - (S1 | T1)
    T2p = T2 - (S1 | T1)
    S2pp, T2pp = S2 - S2p, T2 - T2p
    yield (S1, T1), (S2 | S3, T2 | T3)
    for setI in map(set, combinations(T1 & S2pp, (len(S2pp) - len(T2pp)) // 2)):
        yield (S2pp - setI | S3, T2pp | setI | T3), (S2p | S1 | setI, T2p | T1 - setI)


def relation6BB(S1: set, T1: set, S2: set, T2: set, S3: set, T3: set):
    S2p = S2 - (S1 | T1)
    T2p = T2 - (S1 | T1)
    S2pp, T2pp = S2 - S2p, T2 - T2p
    yield (S2, T2), (S1 | S3, T1 | T3)
    for setI in map(set, combinations(T1 & S2pp, (len(S2pp) - len(T2pp)) // 2)):
        yield (S2pp - setI | S3, T2pp | setI | T3), (S2p | S1 | setI, T2p | T1 - setI)


def H(a_, b_):
    R"""Iterator of (S, T) such that $h_{ST}\in \sH_{ab}$."""
    if b_ - a_ == 0:
        return
    for S, T in Hp(a_ + 1, b_ - 1):
        yield {a_} | S, T | {b_ - 1}


def Hp(a_, b_):
    R"""Iterator of (S, T) such that $h_{ST}\in \sH^\prime_{ab}$."""
    if b_ - a_ == 0:
        yield set(), set()
        return
    for i in range(a_ + 2, b_ + 2, 2):
        for S1, T1 in H(a_, i):
            for S2, T2 in Hp(i, b_):
                yield S1 | S2, T1 | T2


def is_nonzero(S, T):
    return all(i < j for i, j in zip(sorted(S), sorted(T)))


def is_indecomposable(S: set, T: set):
    return is_nonzero(S, T) and is_nonzero(S - {min(S)}, T - {max(T)})


def relations1(j_max):
    for i in range(j_max - 1):
        for j in range(i + 3, j_max + 1):
            yield relation1(i, j)


def relations2(j_max):
    for a1 in range(j_max - 1):
        for b1 in range(a1 + 2, j_max + 1, 2):
            for a2 in range(a1 + 1, b1 + 1, 2):
                for b2 in range(b1 + 1, j_max + 2, 2):
                    for S1, T1 in H(a1, b1):
                        for S2, T2 in H(a2, b2):
                            yield relation2(S1, T1, S2, T2)


def relations3A(j_max):
    for a_ in range(j_max - 2):
        for b_ in range(a_ + 4, j_max + 2, 2):
            N = set(range(a_, b_))
            for j in range(a_ + 1, min(b_ + 1, j_max + 1)):
                for S in map(set, combinations(N, (b_ - a_) // 2 + 1)):
                    T = N - S
                    S1 = {s for s in S if s < j}
                    if S1 and is_indecomposable(S - {max(S1)}, T | {max(S1)}):
                        yield relation3A(j, S, T)


def relations3B(j_max):
    for a_ in range(j_max - 2):
        for b_ in range(a_ + 4, j_max + 2, 2):
            N = set(range(a_, b_))
            for i in range(max(0, a_ - 1), b_ - 1):
                for S in map(set, combinations(N, (b_ - a_) // 2 - 1)):
                    T = N - S
                    T1 = {t for t in T if t > i}
                    if T1 and is_indecomposable(S | {min(T1)}, T - {min(T1)}):
                        yield relation3B(i, S, T)


def relations4A(j_max):
    for a1 in range(j_max - 1):
        for b1 in range(a1 + 4, j_max + 2, 2):
            for a2 in range(a1 + 2, b1, 2):
                for b2 in range(b1, j_max + 2, 2):
                    for S1, T1 in H(a1, b1):
                        for S2, T2 in H(a2, b2):
                            if T2 - T1:
                                yield relation4A(S1, T1, S2, T2)


def relations4B(j_max):
    for a1 in range(j_max - 1):
        for b1 in range(a1 + 2, j_max + 2, 2):
            for a2 in range(a1, b1, 2):
                for b2 in range(b1 + 2, j_max + 2, 2):
                    for S1, T1 in H(a1, b1):
                        for S2, T2 in H(a2, b2):
                            if S1 - S2:
                                yield relation4B(S1, T1, S2, T2)


def relations5(j_max):
    for a1 in range(j_max - 1):
        for b1 in range(a1 + 2, j_max + 2, 2):
            for a2 in range(a1, b1, 2):
                for b2 in range(b1, j_max + 2, 2):
                    for S1, T1 in H(a1, b1):
                        for S2, T2 in H(a2, b2):
                            if b1 - a1 > 2 or b2 - a2 > 2:
                                yield relation5(S1, T1, S2, T2)


def relations6(j_max):
    for a_ in range(j_max - 1):
        for b_ in range(a_ + 4, j_max + 2, 2):
            for (S1, T1), (S2, T2) in combinations(H(a_, b_), 2):
                yield relation6(S1, T1, S2, T2)
            ap, bp = a_ + 1, b_ - 1
            for a1 in range(ap, bp - 5, 2):
                for b1 in range(a1 + 4, bp + 2, 2):
                    for a2, b2 in product(range(a1 + 2, b1, 2), range(b1, bp + 2, 2)):
                        for (S1, T1), (S2, T2) in product(H(a1, b1), H(a2, b2)):
                            if T2 - T1:
                                for (S4, T4), (S5, T5) in product(Hp(ap, a1), Hp(b2, bp)):
                                    for S3, T3 in Hp(a1, a2):
                                        yield sum((h(S6, T6) * h_wrap(S4 | S5 | S7, T4 | T5 | T7)
                                                   for (S6, T6), (S7, T7) in relation6AA(S1, T1, S2, T2, S3, T3)),
                                                  A.zero())
                                    for S3, T3 in Hp(b1, b2):
                                        yield sum((h(S6, T6) * h_wrap(S4 | S5 | S7, T4 | T5 | T7)
                                                   for (S6, T6), (S7, T7) in relation6AB(S1, T1, S2, T2, S3, T3)),
                                                  A.zero())
            for a1 in range(ap, bp - 5, 2):
                for b1 in range(a1 + 2, bp + 2, 2):
                    for a2, b2 in product(range(a1, b1, 2), range(b1 + 2, bp + 2, 2)):
                        for (S1, T1), (S2, T2) in product(H(a1, b1), H(a2, b2)):
                            if S1 - S2:
                                for (S4, T4), (S5, T5) in product(Hp(ap, a1), Hp(b2, bp)):
                                    for S3, T3 in Hp(a1, a2):
                                        yield sum((h(S6, T6) * h_wrap(S4 | S5 | S7, T4 | T5 | T7)
                                                   for (S6, T6), (S7, T7) in relation6BA(S1, T1, S2, T2, S3, T3)),
                                                  A.zero())
                                    for S3, T3 in Hp(b1, b2):
                                        yield sum((h(S6, T6) * h_wrap(S4 | S5 | S7, T4 | T5 | T7)
                                                   for (S6, T6), (S7, T7) in relation6BB(S1, T1, S2, T2, S3, T3)),
                                                  A.zero())


A = GbAlgMod2.load_alg("output/HX9-gens.pickle")
