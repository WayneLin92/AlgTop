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


E2 = GbAlgMod2.load_alg("output/HX7.pickle")
A = GbAlgMod2.new_alg(key=key)
A.add_gens(zip(E2.gen_names))
