from algebras.operations import Steenrod
from algebras.groebner import FreeModuleMod2 as FM
from algebras import linalg
from spec.specseq import SpecSeq
from typing import List, Tuple


class Ext_A:
    def __init__(self, s_max, t_max):
        # todo: abandon s_max
        self.s_max = s_max
        self.t_max = t_max
        self.h = [[] for _ in range(s_max + 1)]  # type: List[List[Tuple[int, FM]]]

    def compute_minimal(self):
        """ return Ext_R(k, k) """
        FM.set_ring(Steenrod)

        s_max, t_max = self.s_max, self.t_max
        self.h[0].append((0, FM.gen(1)))  # (internal_deg, diff)

        # initialize the augmented ideal of ring
        ideal = [None]  # type: List[Tuple[Steenrod]]
        for t in range(1, t_max + 1):
            ideal.append(tuple(Steenrod.basis(t)))

        # compute the resolution by internal_deg
        for t in range(1, t_max + 1):
            kernel_t = linalg.VectorSpaceMod2(r * FM.gen("a_{0,0}") for r in ideal[t])
            for s in range(1, s_max + 1):
                print("(s, t)=({}, {})".format(s, t))
                my_map = linalg.LinearMapKMod2()
                for index_gen in range(len(self.h[s])):
                    t_gen, gen_kernel = self.h[s][index_gen]
                    my_map.add_maps((sq * FM.gen("a_{{{}, {}}}".format(s, index_gen)), sq * gen_kernel)
                                    for sq in ideal[t-t_gen])
                for gen in (kernel_t / my_map.image).basis(FM):
                    self.h[s].append((t, gen))
                kernel_t = my_map.kernel

    def __str__(self):
        result = ""
        for s in range(len(self.h)):
            for j in range(len(self.h[s])):
                result += "$a_{{{}, {}}} ({}, {}) \\to {}$\\\\\n".\
                    format(s, j, self.h[s][j][0] - s, s, self.h[s][j][1])
            result += "\\\\\n"
        return result

    def get_spec(self, x_max, y_max) -> SpecSeq:
        spec = SpecSeq(x_max, y_max, "Adams Discrete")
        for s in range(len(self.h)):
            for j in range(len(self.h[s])):
                t = self.h[s][j][0]
                if s <= y_max and t - s <= x_max:
                    spec.add_single_gen("a_{{{}, {}}}".format(s, j), (s, t))
        return spec

    # TODO: save and load


# tests -----------------------------
def test_ext():
    ext = Ext_A(15, 30)
    ext.compute_minimal()
    spec = ext.get_spec(30, 15)
    print(ext)
    spec.draw()


def test():
    FM.set_ring(Steenrod)
    print(Steenrod.gen(2) * FM.gen("x"))


if __name__ == "__main__":
    # test_amod()
    test_ext()
    # test()
