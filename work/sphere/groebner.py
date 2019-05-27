import sympy as sp


def alg_B(n_max):
    B = {}
    for i in range(n_max):
        for j in range(i + 1, n_max + 1):
            B[(i, j)] = sp.Symbol(f'B^{i}_{j}')
    rels = []
    for d in range(2, n_max + 1):
        for i in range(n_max + 1 - d):
            j = i + d
            rels.append(sum((B[(i, k)] * B[(k, j)] for k in range(i + 1, j))))
    fns = sp.groebner(rels).args[0]
    result = []
    for f in fns:
        result.append(f"${sp.latex(f.as_expr())}$\\\\")
    result.sort(key=len)
    for s in result:
        if "0" in s:
            print(s)


if __name__ == "__main__":
    alg_B(7)
