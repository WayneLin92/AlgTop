"""Class of spectral sequences."""
# TODO: log file

from GUI.draw import draw_ss

import copy
import collections
import itertools
from typing import Union, Optional, List, Dict, Tuple, Set, Callable
from algebras import linalg, mymath
from algebras.constructions import AugAlgMod2

TYPE_IMAGE = 0
TYPE_KERNEL = 1
TYPE_DIFF = 2
TYPE_TBD = 3


class MyPoly:
    pass


class SpecSeq:
    """Spectral sequence.

    All the internal degree has data type mymath.Deg.
    """
    MyTuple = collections.namedtuple('MyTuple', ('p_max', 'q_max', 'Alg', 'basis', 'diff_map'))

    def __init__(self, p_max: int, q_max: int, *, starting_page=2, func=None, Alg=None):
        """func is a function indicating the direction of the differential for each page."""
        self.data = [self.MyTuple(p_max, q_max, Alg or AugAlgMod2.new_alg(), None, linalg.GradedLinearMapMod2())]
        self._locked = False  # flag for using Alg.add_gen, Alg.add_map
        self._starting_page = starting_page
        self._func_deg_diff = func if func else lambda r: (r, r + 1)  # type: Callable[[int], Tuple[int, int]]

    @property
    def page(self):
        """Get the page number of the last page."""
        return self._starting_page + len(self.data) - 1

    @property
    def deg_diff(self):
        """Get the direction of the differential of the last page."""
        return mymath.Deg(self._func_deg_diff(self.page))

    def present(self, deg):
        pass

    def draw(self):
        draw_ss(self)

    def get_basis(self):
        """Get the basis of self.data[-1].Alg. Return type: Iterator[deg, Iterator[tuple]]."""
        p_max, q_max = self.data[-1].p_max, self.data[-1].q_max
        d_max = mymath.Deg((p_max, q_max))
        R = self.data[-1].Alg
        basis_mon = sorted(R.basis_mons_max(d_max), key=R.deg_mon)
        return itertools.groupby(basis_mon, key=R.deg_mon)

    # GUI interface
    def get_bullets_and_arrows(self) -> Tuple[Dict[tuple, list], List[Tuple[tuple, int]]]:
        """Get a representation of self for drawing the spectral sequence."""
        R = self.data[-1].Alg
        linmap = self.data[-1].diff_map
        if self._locked:
            basis = self.data[-1].basis  # type: Dict[mymath.Deg, Set[tuple]]
            bullets = {}
            arrows = []
            for deg in basis:
                deg_src = deg - self.deg_diff
                vs_image = linmap.image(deg_src)
                vs_kernel = linmap.kernel(deg)
                list_inv_image = linmap.data[deg].g
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
            basis = self.get_basis()
            bullets = {d: [(R(m), TYPE_TBD) for m in v] for d, v in basis}
            return bullets, []

    def add_gen(self, k, deg: tuple):
        if self._locked:
            print("Error: the algebra is locked.")
            return
        deg = mymath.Deg(deg)
        self.data[0].Alg.add_gen(k, deg)

    def add_rel(self, rel):
        if self._locked:
            print("Error: the algebra is locked.")
            return
        self.data[0].Alg.add_rel(rel)

    def add_diff(self, src, tgt):
        deg_diff = tgt.deg() - src.deg()
        if deg_diff != self.deg_diff:
            print("Error: wrong degree of differential.")
            return
        self.data[-1].diff_map.add_map(src, tgt)

        if not self._locked:
            self._locked = True
            self.data[-1].basis = {}
            pass

    def add_diff1(self, poly_saddr1, poly_saddr2, b_log=True):
        # type: (Union[MyPoly, tuple, None], Union[MyPoly, tuple, None], bool) -> Optional[MyPoly]
        """
        partly interface for draw.py
        return the source poly of the differential
        """
        if type(poly_saddr1) is tuple:  # interface for draw.py
            surf_deg1, index1 = poly_saddr1
            deg1 = self.sd2deg(surf_deg1)
            index1 += self.indices[-1][deg1][IINDEX_IMAGE]
            bul1 = self.entries[-1][deg1][index1]
            poly_saddr1 = bul1['poly']
            if poly_saddr2 is not None:
                surf_deg2, index2 = poly_saddr2
                deg2 = self.sd2deg(surf_deg2)
                index2 += self.indices[-1][deg2][IINDEX_IMAGE]
                poly_saddr2 = self.entries[-1][deg2][index2]['poly']
            else:
                deg2 = None
                poly_saddr2 = MyPoly.zero()
        else:
            deg1 = poly_saddr1.deg()
            deg2 = poly_saddr2.deg()
            bul1 = self._get_bullet(poly_saddr1)

        if bul1 is not None and bul1['diff'] is not None:
            if bul1['diff'] != poly_saddr2:
                raise SSBulletsError('Contradictory differentials: d({})={} vs d({})={}'.format(
                    poly_saddr1, poly_saddr2, poly_saddr1, bul1['diff']))
            else:
                if b_log:
                    print('{}: repeating differentials!'.format(deg1))
                return None
        if deg2 is None or deg2 == self._diff_target(deg1):
            bullet1 = {'poly': poly_saddr1, 'base': None, 'mon': None, 'diff': poly_saddr2, 'type': TYPE_DIFF}
            msg_index1 = self.add_bullet(bullet1, deg1)
            if b_log:
                nb_msg(msg_index1, deg1)
            if deg2 is not None:
                bullet2 = {'poly': poly_saddr2, 'base': None, 'mon': None, 'diff': MyPoly.zero(),
                           'inv_diff': poly_saddr1, 'base_inv_diff': None, 'type': TYPE_IMAGE}
                self.add_bullet(bullet2, deg2)
            return bullet1['poly'] if msg_index1 == NB_LESS_TBD else None
        return None

    def deduce_init(self):
        if not self.multiplicative:
            return
        for i, j in self._get_degs("trivial_diff"):
            for bullet in self._get_nonrel_bullets((i, j)):
                self.add_diff(bullet['poly'], MyPoly.zero(), False)
        for i, j in self._get_degs("square"):
            for bullet in self._get_nonrel_bullets((i, j)):
                self.add_diff(bullet['poly'].frob(), MyPoly.zero(), False)

    def deduce_new_diff(self, poly):
        if not self.multiplicative:
            return
        deg, index = self._get_addr(poly)
        index += self.indices[-1][deg][IINDEX_IMAGE]
        bullet = self.entries[-1][deg][index]
        if bullet['type'] != TYPE_DIFF:
            raise SSValueError("{}: bullet should be of type diff, got ({}, type{}) instead".
                               format(deg, bullet['poly'], bullet['type']))
        for i, j in self._get_degs("src_diff", deg):
            index_image = self.indices[-1][(i, j)][IINDEX_IMAGE]
            index_tbd = self.indices[-1][(i, j)][IINDEX_TBD]
            for k in range(index_image, index_tbd):
                bul = self.entries[-1][(i, j)][k]
                self.add_diff(bullet['poly'] * bul['poly'],
                              bullet['diff'] * bul['poly'] + bullet['poly'] * bul['diff'], False)

    # TODO def deduce_no_choice
    # TODO add relations
    def homology(self, deg):
        result = copy.deepcopy(self.entries[-1][deg])
        for bullet in result:
            if bullet['type'] == TYPE_REL:
                pass
            elif bullet['type'] == TYPE_IMAGE:
                del bullet['inv_diff']
                del bullet['base_inv_diff']
                bullet['type'] = TYPE_REL
            elif bullet['type'] == TYPE_DIFF:
                if not bullet['diff']:
                    bullet['diff'] = None
                    bullet['base'] = None
                    del bullet['base_diff']
                    bullet['type'] = TYPE_TBD
                else:
                    bullet['type'] = None
            elif bullet['type'] == TYPE_TBD:
                bullet['base'] = None
        result = [bullet for bullet in result if bullet['type'] is not None]
        return result

    @staticmethod
    def gen_indices(bullets):
        result = []
        i = 0
        while i < len(bullets) and bullets[i]['type'] == TYPE_REL:
            i += 1
        result.append(i)
        while i < len(bullets) and bullets[i]['type'] == TYPE_IMAGE:
            i += 1
        result.append(i)
        while i < len(bullets) and bullets[i]['type'] == TYPE_DIFF:
            i += 1
        result.append(i)
        return result

    def new_page(self):
        self.p_max.append(self.p_max[-1] - self.page)
        self.q_max.append(self.q_max[-1] - self.page + 1)
        entries_page = dict((deg, self.homology(deg)) for deg in self.entries[-1]
                            if deg[0] <= self.p_max[-1] and deg[1] <= self.q_max[-1])
        indices_page = dict((deg, self.gen_indices(bullets)) for deg, bullets in entries_page.items())
        self.entries.append(entries_page)
        self.indices.append(indices_page)
        self.page += 1
        self._gen_basis_new_page()
        self.deduce_init()

# 623
