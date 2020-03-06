"""Class of spectral sequences."""
# TODO: next_page
# TODO: improve GUI
import collections
import itertools
import operator
from typing import List, Dict, Tuple, Set, Callable
from algebras import linalg, mymath
from algebras.groebner import GbAlgMod2

TYPE_IMAGE = 0
TYPE_KERNEL = 1
TYPE_DIFF = 2
TYPE_TBD = 3


class SpecSeq:
    """Spectral sequence.

    All the internal degrees have data type mymath.Deg.
    """
    MyTuple = collections.namedtuple('MyTuple', ('mask', 'Alg', 'basis', 'diff_map'))

    def __init__(self, p_max: int, q_max: int, *, starting_page=2, func=None, Alg=None):
        """func is a function indicating the direction of the differential for each page."""
        mask = set(map(mymath.Deg, itertools.product(range(p_max + 1), range(q_max + 1))))
        self.data = [self.MyTuple(mask, Alg or GbAlgMod2.new_alg(unit_deg=mymath.Deg((0, 0))),
                                  {}, linalg.GradedLinearMapMod2())]
        self._d_max = mymath.Deg((p_max, q_max))
        self._initialized = False  # flag for using Alg.add_gen, Alg.add_map
        self._starting_page = starting_page
        self._func_deg_diff = func if func else lambda r: (r, 1 - r)  # type: Callable[[int], Tuple[int, int]]

    # getters -----------------
    @property
    def page(self):
        """Get the page number of the last page."""
        return self._starting_page + len(self.data) - 1

    @property
    def deg_diff(self):
        """Get the direction of the differential of the last page."""
        return mymath.Deg(self._func_deg_diff(self.page))

    def get_basis(self):
        """Get the basis of self.data[-1].Alg."""
        R, mask = self.data[-1].Alg, self.data[-1].mask
        result = {}
        for d, m in R.basis_mons_max(self._d_max):
            if d in mask:
                if d in result:
                    result[d].append(m)
                else:
                    result[d] = [m]
        return result.items()

    def present(self, deg):
        pass

    # setters -------------------
    def add_gen(self, k, deg: tuple):
        if self._initialized:
            print("Error: the algebra can no longer be changed.")
            return
        deg = mymath.Deg(deg)
        return self.data[0].Alg.add_gen(k, deg)

    def add_rel(self, rel):
        if self._initialized:
            print("Error: the algebra can no longer be changed.")
            return
        self.data[0].Alg.add_rel(rel)

    def _add_diff(self, src, tgt):
        """Add a single differential."""
        deg1, deg2 = src.deg(), tgt.deg()
        mask = self.data[-1].mask
        if deg1 is None:
            if deg2 is not None:
                raise ValueError("incompatible differential")
            return
        if deg2 is None:
            deg2 = deg1 + self.deg_diff
            if deg1 in mask and (deg2 in mask or deg2[0] < 0 or deg2[1] < 0):
                self.data[-1].diff_map.add_map(src, tgt)
            return
        deg_diff = deg2 - deg1
        if deg_diff != self.deg_diff:
            print(f"Error: wrong degree of differential. Expected: {self.deg_diff}, got: {deg_diff}.")
            return
        if deg1 in mask and deg2 in mask:
            self.data[-1].diff_map.add_map(src, tgt)

    def add_diff(self, src, tgt):
        if self._initialized:
            R, linmap = self.data[-1].Alg, self.data[-1].diff_map
            domain = [x for d in linmap.data for x in linmap.domain(d).basis(R)]
            for x in domain:
                src1 = src * x
                tgt1 = tgt * x + src * linmap.f(x)
                self._add_diff(src1, tgt1)
                self._add_diff(tgt1, tgt1.zero())

    def init_diff(self):
        """Add differential for degree reason and Frobenius."""
        if not self._initialized:
            self._initialized = True
            self.data[0].basis.update((d, set(g)) for d, g in self.get_basis())
            mask, R, basis = self.data[0].mask, self.data[0].Alg, self.data[0].basis
            squares = set()
            for d in mask:
                if d * 2 in mask and d in basis:
                    for m in basis[d]:
                        self._add_diff(R(m).square(), R.zero())
                        squares.add((i * 2 for i in m))
            for d in mask:
                deg2 = d + self.deg_diff
                if deg2[0] < 0 or deg2[1] < 0 and d in basis:
                    for m in basis[d]:
                        self.add_diff(R(m), R.zero())

    # todo: def deduce_no_choice
    def homology(self, deg):
        pass

    def new_page(self):
        pass

    # GUI interface -------------------
    def get_bullets_and_arrows(self) -> Tuple[Dict[tuple, list], List[Tuple[tuple, int]]]:
        """Get a representation of self for drawing the spectral sequence."""
        R = self.data[-1].Alg
        linmap = self.data[-1].diff_map
        if self._initialized:
            basis = self.data[-1].basis  # type: Dict[mymath.Deg, Set[tuple]]
            bullets = {}
            arrows = []
            for deg in basis:
                deg_src = deg - self.deg_diff
                vs_image = linmap.image(deg_src)
                vs_kernel = linmap.kernel(deg)
                list_inv_image = linmap.data[deg].g if deg in linmap.data else []
                set_TBD = basis[deg] - set(linmap.domain(deg).get_mons())

                bullets[deg] = [(x, TYPE_IMAGE) for x in vs_image.basis(R)]
                bullets[deg].extend((x, TYPE_KERNEL) for x in (vs_kernel / vs_image).basis(R))
                bullets[deg].extend((x, TYPE_DIFF) for x in map(R, list_inv_image))
                bullets[deg].extend((x, TYPE_TBD) for x in map(R, sorted(set_TBD)))
                dim_kernel = vs_kernel.dim
                num_arrow = len(list_inv_image)
                arrows.extend(((deg, dim_kernel + i), (deg + self.deg_diff, i)) for i in range(num_arrow))
            return bullets, arrows
        else:
            bullets = {d: [(R(m), TYPE_TBD) for m in sorted(g)] for d, g in self.get_basis()}
            return bullets, []

    def draw(self):
        import GUI.draw
        GUI.draw.draw_ss(self)
        del GUI


def test():
    spec = SpecSeq(10, 10, starting_page=1, func=lambda r: (r, -r))
    h0 = spec.add_gen("h0", (1, 0))
    h1 = spec.add_gen("h1", (1, 0))
    b02 = spec.add_gen("b02", (2, 0))
    h2 = spec.add_gen("h2", (0, 1))
    R12 = spec.add_gen("R12", (0, 2))
    R03 = spec.add_gen("R03", (0, 3))
    spec.add_rel(h0 * h1)
    spec.init_diff()
    # spec.add_diff(h2, h2.zero())
    # spec.add_diff(R12, h1*h2)
    # spec.add_diff(R03, h0*R12)
    spec.draw()
    return spec


if __name__ == "__main__":
    test()

# 623, 193, 184
