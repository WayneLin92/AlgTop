R"""Check relations in HX_{n, m}."""
from algebras.groebner import GbAlgMod2
from typing import Union
import itertools


E1 = GbAlgMod2.new_alg(key=lambda _m: [-i for i in _m])


def R(S: Union[int, tuple], T: Union[int, tuple]):
    if type(S) is int:
        return E1.gen(f"R_{{{S}{T}}}")
    assert len(S) == len(T)
    S, T = sorted(S), sorted(T)
    result = E1.zero()
    for T1 in itertools.permutations(T):
        if all(t - s > 0 for s, t in zip(S, T1)):
            pro = E1.unit()
            for x in map(R, S, T1):
                pro *= x
            result += pro
    return result


def b(S, T):
    return R(S, T) ** 2


def h(S: Union[int, tuple], T: tuple = None):
    """Return h_i(S)."""
    if type(S) is int:
        if T is None:
            return R(S, S + 1)
        i, S1 = S, T
        k = len(S1)
        seq = set(range(i, i + 2 * k + 2))
        S = {i + s for s in S1} | {i}
        T = seq - S
        assert len(S) + len(T) == 2 * k + 2
    S, T = sorted(S), sorted(T)
    result = E1.zero()
    for T1 in itertools.permutations(T):
        if all(t - s > 0 for s, t in zip(S, T1)):
            prod = E1.unit()
            for s, t in zip(S, T1):
                prod *= R(s, t)
            result += prod
    return result


def gen_MayE1(m, n_max):
    """Return $X_{n_max, m}/im d$."""
    gens = []
    for i in range(n_max):
        for j in range(i + 1, n_max + 1):
            if i > 0 or j <= m:
                gens.append((f"R_{{{i}{j}}}", 2 ** j - 2 ** i, j - i))
    gens.sort(key=lambda _x: _x[2])

    E1.add_gens(gens)

    rels = []
    for d in range(2, n_max + 1):
        for i in range(n_max + 1 - d):
            j = i + d
            if i > 0 or j <= m:
                rel = sum((R(i, k) * R(k, j) for k in range(i + 1, j)), E1.zero())
                rels.append(rel)
    E1.add_rels(rels)
