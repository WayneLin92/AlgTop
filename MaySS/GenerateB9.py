"""Generate B9."""
from algebras.groebner import GbAlgMod2


def key(mon):
    return [-_i for _i in mon]


E1 = GbAlgMod2.new_alg(key=key)


def R(S, T):
    return E1.gen(f"R_{{{S}{T}}}")


if __name__ == "__main__":
    n_max = 9
    gens = []
    for i in range(n_max):
        for j in range(i + 1, n_max + 1):
            gens.append((f"R_{{{i}{j}}}", 2 ** j - 2 ** i, j - i))
    gens.sort(key=lambda _x: _x[2])
    E1.add_gens(gens)

    rels = []
    for d in range(2, n_max + 1):
        for i in range(n_max + 1 - d):
            j = i + d
            print(i, j)
            rel = sum((R(i, k) * R(k, j) for k in range(i + 1, j)), E1.zero())
            rels.append(rel)
            E1.add_rel(rel)

    E1.save_alg('../output/B9.pickle')
