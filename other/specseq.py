import algebras.polynomials
import copy
from numpy import array
from GUI.draw import draw_ss
from GUI.constants import *
from typing import Union, Optional, Iterator, Tuple, List, Dict

TYPE_REL = 0
TYPE_IMAGE = 1
TYPE_DIFF = 2
TYPE_TBD = 3

IINDEX_IMAGE = 0
IINDEX_DIFF = 1
IINDEX_TBD = 2
IINDEX_NUM = 3
IINDEX_INSERT = {TYPE_REL: IINDEX_IMAGE, TYPE_IMAGE: IINDEX_DIFF, TYPE_DIFF: IINDEX_TBD, TYPE_TBD: None}

NB_REDUNDANT_REL = 0
NB_REDUNDANT_IMAGE = 1
NB_REDUNDANT_DIFF = 2
NB_REDUNDANT_IMAGE1 = 3
NB_REDUNDANT_DIFF1 = 4
NB_LESS_TBD = 5
NB_NEW = 6

SS_TYPE_SERRE_COHOMOLOGY = 0
SS_TYPE_ADAMS = 1
SS_TYPE_DICT = {"Serre-Cohomology": SS_TYPE_SERRE_COHOMOLOGY,
                "Adams": SS_TYPE_ADAMS}


def nb_msg(msg_index, deg):
    if msg_index == NB_REDUNDANT_REL:
        print("{}: Extra relation!".format(deg))
    elif msg_index == NB_REDUNDANT_IMAGE:
        print("{}: There is a relation between images".format(deg))
    elif msg_index == NB_REDUNDANT_DIFF:
        print("{}: Redundant differential".format(deg))
    elif msg_index == NB_REDUNDANT_IMAGE1:
        print("{}: 2. There is a relation between images".format(deg))
    elif msg_index == NB_REDUNDANT_DIFF1:
        print("{}: 2. Redundant differential!".format(deg))
    elif msg_index == NB_LESS_TBD:
        pass
        # print("{}: 2. One less TBD!".format(deg))


# TODO: log file


class SSError(Exception):
    pass


class SSBulletsError(SSError):
    pass


class SSValueError(SSError):
    pass


class SSClassError(SSError):
    pass


def poly_proj_mod(poly, bullets) -> Tuple[list, 'MyPoly']:
    """ return (proj, res) where poly = proj + res """
    proj = []
    res = poly.copy()
    for bullet in bullets:
        if bullet['mon'] in res.data:
            coeff = res.data[bullet['mon']]
            proj.append(coeff)
            res -= bullet['base'] * coeff
        else:
            proj.append(0)
    return proj, res


def poly_mod(poly, bullets):
    """ return poly - proj(poly, bullets) """
    result = poly.copy()
    for bullet in bullets:
        if bullet['mon'] in result.data:
            result -= bullet['base'] * result.data[bullet['mon']]
    return result


class MyPoly(algebras.polynomials.PolyAnyVarModP):
    # -- PolyModP ---------
    @classmethod
    def gen(cls, key, deg: tuple = (1, 0)):
        if key in cls.dict_deg_gen and cls.dict_deg_gen[key] != deg:
            raise SSValueError("Warning: the degree of {} is changed from {} to {}".
                               format(key, cls.dict_deg_gen[key], deg))
        cls.dict_deg_gen[key] = deg
        return cls(((key, 1),))

    def deg_mon(self, mon) -> tuple:
        try:
            deg = sum((self.deg_gen(gen) * exp for gen, exp in mon), array([0, 0]))
        except ValueError:
            print(mon)
            raise ValueError
        return deg[0], deg[1]

    def deg_gen(self, key) -> array:
        return array(self.dict_deg_gen[key])

    # methods
    def get_mon(self):  # ######
        """ return a monomial """
        if len(self.data) > 0:
            return max(self.data)
        else:
            return None

    def get_item(self):  # ######
        """ return a monomial with its coefficient """
        if len(self.data) > 0:
            return max(self.data.items())
        else:
            return None, None

    @staticmethod
    def prod_mons(mon1, mon2):
        """ return the product of two monomials as a monomial """
        mon_pro = dict(mon1)
        for gen, exp in mon2:
            if gen in mon_pro:
                mon_pro[gen] += exp
            else:
                mon_pro[gen] = exp
        return tuple(sorted(mon_pro.items()))


class SpecSeq:
    """
    self.gen_deg[id] = (deg_x, deg_y)
    self.entries[page-2][(x,y)] is a list of bullets {'poly': ?, 'base': ?, 'mon': ?, 'diff': ?, 'type': ?}
    self.indices[page-2][(x,y)] is [index_image, index_diff, index_TBD]
    """

    def __init__(self, x_max: int, y_max: int, ss_type: str = "Serre-Cohomology"):
        self.x_max = [x_max]
        self.y_max = [y_max]
        if ss_type == "Serre-Cohomology":
            self.entries = [{(0, 0): [{'poly': MyPoly.unit(), 'base': MyPoly.unit(), 'mon': (),
                                       'diff': MyPoly.zero(), 'type': TYPE_DIFF, 'base_diff': MyPoly.zero()}]}]
        else:
            self.entries = [{}]  # type: List[Dict[tuple, dict]]
        self.indices = [{(0, 0): [0, 0, 1]}]
        self.page = 2
        strings = ss_type.split()
        self.ss_type = SS_TYPE_DICT[strings[0]]
        self.multiplicative = False if len(strings) > 1 and strings[1] == "Discrete" else True

    def present(self, keys: Iterator[str], deg: Optional[tuple] = None):
        if deg is None:
            for deg, bullets in sorted(self.entries[self.page-2].items()):
                # noinspection PyTypeChecker
                print("{}: {}".format(deg, [[str(bullet[key]) if key is not 'mon' else str(MyPoly({bullet[key]: 1}))
                                             for key in keys] for bullet in bullets]))
        else:
            bullets = self.entries[self.page-2][deg]
            print("{}: {}".format(deg, [[str(bullet[key]) if key is not 'mon' else str(MyPoly({bullet[key]: 1}))
                                         for key in keys] for bullet in bullets]))

    def draw(self):
        draw_ss(self)

    # interface functions for draw.py
    def sd2deg(self, surf_deg: Optional[tuple]) -> Optional[tuple]:
        """ interface for draw.py """
        if surf_deg is None:
            return None
        if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
            return surf_deg
        elif self.ss_type == SS_TYPE_ADAMS:
            return surf_deg[1], surf_deg[0] + surf_deg[1]
        else:
            raise SSClassError

    def deg2sd(self, deg: Optional[tuple]) -> Optional[tuple]:
        """ interface for draw.py """
        if deg is None:
            return None
        if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
            return deg
        elif self.ss_type == SS_TYPE_ADAMS:
            return deg[1] - deg[0], deg[0]
        else:
            raise SSClassError

    def get_surf_degs(self) -> Iterator:
        for deg in self.entries[self.page - 2]:
            yield self.deg2sd(deg)

    def get_surf_addr(self, poly, surf_deg=None) -> Optional[tuple]:
        """ interface for draw.py """
        addr = self._get_addr(poly, surf_deg)
        if addr is None:
            return None
        else:
            deg, index = addr
            return self.deg2sd(deg), index

    def get_poly(self, surf_addr):
        """ interface for draw.py """
        surf_deg, index = surf_addr
        deg = self.sd2deg(surf_deg)
        index_image = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        return self.entries[self.page - 2][deg][index_image + index]['poly']

    def get_bullet_color(self, surf_addr):
        """ interface for draw.py """
        surf_deg, index = surf_addr
        deg = self.sd2deg(surf_deg)
        index_image = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        bullet = self.entries[self.page - 2][deg][index_image + index]
        return COLOR_BULLET_TBD if bullet['type'] == TYPE_TBD else COLOR_BULLET_NON_TBD

    def get_num_nonrel_bullets(self, surf_deg) -> int:
        deg = self.sd2deg(surf_deg)
        num_rel = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        num_total = len(self.entries[self.page - 2][deg])
        return num_total - num_rel

    def get_arrows(self, expansion):
        """
        interface for draw.py
        Iterator of arrows
        """
        for deg, bullets in self.entries[self.page - 2].items():
            index_image = self.indices[self.page - 2][deg][IINDEX_IMAGE]
            index_diff = self.indices[self.page - 2][deg][IINDEX_DIFF]
            index_tbd = self.indices[self.page - 2][deg][IINDEX_TBD]
            for i in range(index_diff, index_tbd):
                diff = bullets[i]['diff']
                if diff:
                    yield (self.deg2sd(deg), i-index_image), self.get_surf_addr(diff)
                    if self.get_num_nonrel_bullets(self.deg2sd(deg)) > 9 and \
                            self.get_num_nonrel_bullets(self.deg2sd(self._diff_target(deg))) > 9 \
                            and self.deg2sd(deg) not in expansion and \
                            self.deg2sd((deg[0] + self.page, deg[1] - self.page + 1)) not in expansion:
                        break
    # end of interface

    def _get_degs(self, type_bullets: str, truncate_deg: tuple = (0, 0)) -> Iterator:
        """ Iterator of degrees """
        if type_bullets == "trivial_diff":
            if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
                for i in range(0, self.x_max[-1] + 1):
                    for j in range(0, min(self.page - 1, self.y_max[-1] + 1)):
                        if (i, j) in self.entries[-1]:
                            yield i, j
            elif self.ss_type == SS_TYPE_ADAMS:
                for s in range(0, self.y_max[-1] + 1):
                    if (s, s) in self.entries[-1]:
                        yield s, s
        elif type_bullets == "square":
            p = MyPoly.get_prime()
            if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
                for i in range(0, self.x_max[-1] // p + 1):
                    for j in range((self.page - 1 + p - 1) // p, self.y_max[-1] // p + 1):
                        if (i, j) in self.entries[-1]:
                            yield i, j
            elif self.ss_type == SS_TYPE_ADAMS:
                for s in range(0, self.y_max[-1] // p + 1):
                    for t in range(s + 1, s + self.x_max[-1] // p + 1):
                        if (s, t) in self.entries[-1]:
                            yield s, t
        if type_bullets == "src_diff":
            if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
                for i in range(0, self.x_max[-1] - truncate_deg[0] - self.page + 1):
                    for j in range(0, self.y_max[-1] - truncate_deg[1] + 1):
                        if (i, j) in self.entries[-1]:
                            yield i, j
            elif self.ss_type == SS_TYPE_ADAMS:
                for s in range(0, self.y_max[-1] - truncate_deg[0] - self.page + 1):
                    for t in range(s, s + self.x_max[-1] - (truncate_deg[1] - truncate_deg[0]) + 1):
                        if (s, t) in self.entries[-1]:
                            yield s, t
        if type_bullets == "bullet":
            if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
                for i in range(0, self.x_max[-1] - truncate_deg[0] + 1):
                    for j in range(0, self.y_max[-1] - truncate_deg[1] + 1):
                        if (i, j) in self.entries[-1]:
                            yield i, j
            elif self.ss_type == SS_TYPE_ADAMS:
                for s in range(0, self.y_max[-1] - truncate_deg[0] + 1):
                    for t in range(s, s + self.x_max[-1] - (truncate_deg[1] - truncate_deg[0]) + 1):
                        if (s, t) in self.entries[-1]:
                            yield s, t

    def _inside_range(self, deg):
        """ Test if deg is in the region determined by x_max, y_max """
        if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
            return 0 <= deg[0] <= self.x_max[-1] and 0 <= deg[1] <= self.y_max[-1]
        elif self.ss_type == SS_TYPE_ADAMS:
            return 0 <= deg[0] <= self.y_max[-1] and 0 <= deg[1] - deg[0] <= self.x_max[-1]

    def _get_addr(self, poly, surf_deg=None) -> Optional[tuple]:
        """ interface for draw.py """
        deg = self.sd2deg(surf_deg)
        if deg is None:
            deg = poly.deg()
        bullets = self.entries[self.page - 2][deg]
        index_image = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        num_total = len(bullets)
        for i in range(index_image, num_total):
            if bullets[i]['poly'] == poly:
                return deg, i - index_image
        return None

    def _get_bullet(self, poly: MyPoly, deg=None):
        """ modified from get_surf_addr """
        if deg is None:
            deg = poly.deg()
        bullets = self.entries[self.page - 2][deg]
        index_image = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        num_total = len(bullets)
        for i in range(index_image, num_total):
            if bullets[i]['poly'] == poly:
                return bullets[i]
        return None

    def _get_nonrel_bullets(self, deg):
        """ Iterator for nonrel bullets """
        num_rel = self.indices[self.page - 2][deg][IINDEX_IMAGE]
        num_total = len(self.entries[self.page - 2][deg])
        for i in range(num_rel, num_total):
            yield self.entries[self.page - 2][deg][i]

    def _diff_target(self, deg_src):
        if self.ss_type == SS_TYPE_SERRE_COHOMOLOGY:
            return deg_src[0] + self.page, deg_src[1] - self.page + 1
        elif self.ss_type == SS_TYPE_ADAMS:
            return deg_src[0] + self.page, deg_src[1] + self.page - 1

    def _gen_basis_new_bullet(self, deg: tuple, index_nb: int) -> int:
        """ executed when a new bullet is added """
        result = None
        bullets = self.entries[-1][deg]
        indices = self.indices[-1][deg]

        i = index_nb
        # deal with the new added bullet
        proj, poly_res = poly_proj_mod(bullets[i]['poly'], bullets[:i])
        mon, coeff = poly_res.get_item()
        if mon is None:
            if bullets[i]['type'] == TYPE_REL:
                result = NB_REDUNDANT_REL
            elif bullets[i]['type'] == TYPE_IMAGE:  #
                result = NB_REDUNDANT_IMAGE
                bul_inv = self._get_bullet(bullets[i]['inv_diff'])
                if bul_inv is not None:
                    bul_inv['poly'] -= sum(bullets[j]['base_inv_diff'] * proj[j]
                                           for j in range(indices[IINDEX_IMAGE], i))
                    bul_inv['diff'] = MyPoly.zero()
            elif bullets[i]['type'] == TYPE_DIFF:
                diff = sum(bullets[j]['base_diff'] * proj[j] for j in range(indices[IINDEX_DIFF], i))
                if diff != bullets[i]['diff']:
                    raise SSBulletsError('{}: Contradictory differentials: d({})={} vs d({})={}'.
                                         format(deg, bullets[i]['poly'], diff, bullets[i]['poly'], bullets[i]['diff']))
                else:
                    result = NB_REDUNDANT_DIFF
            else:
                raise SSBulletsError("{}:  Unexpected new tbd bullet!".format(deg))
            del bullets[i]
            for j in range(len(indices)):
                if indices[j] > i:
                    indices[j] -= 1
            return result
        else:
            poly = poly_res * MyPoly.INV[coeff]
            bullets[i]['base'] = poly
            bullets[i]['mon'] = mon
            if bullets[i]['type'] == TYPE_DIFF:
                bullets[i]['base_diff'] = (bullets[i]['diff'] -
                                           sum(bullets[j]['base_diff'] * proj[j]
                                               for j in range(indices[IINDEX_DIFF], i))) * MyPoly.INV[coeff]
            elif bullets[i]['type'] == TYPE_IMAGE:
                bullets[i]['base_inv_diff'] = (bullets[i]['inv_diff'] -
                                               sum(bullets[j]['base_inv_diff'] * proj[j]
                                                   for j in range(indices[IINDEX_IMAGE], i))) * MyPoly.INV[coeff]

        # deal with the bullets after the new bullet
        index_nb = i  # index of the new bullet
        i = index_nb + 1
        while i < len(bullets) and bullets[index_nb]['mon'] not in bullets[i]['base'].data:
            i += 1
        if i >= len(bullets) and index_nb < len(bullets) - 1:
            raise SSBulletsError("{}: Expecting a relation".format(deg))
        while i < len(bullets):
            proj, poly_res = poly_proj_mod(bullets[i]['base'], bullets[index_nb:i])  #
            mon, coeff = poly_res.get_item()
            if mon is None:
                # not possible to be TYPE_REL
                if bullets[i]['type'] == TYPE_IMAGE:  # this case is very unusual
                    result = NB_REDUNDANT_IMAGE1
                elif bullets[i]['type'] == TYPE_DIFF:
                    if bullets[i]['base_diff']:
                        raise SSBulletsError("{}: 2. Contradictory differentials!".format(deg))
                    else:
                        result = NB_REDUNDANT_DIFF1
                else:
                    result = NB_LESS_TBD
                del bullets[i]
                for j in range(len(indices)):
                    if indices[j] > i:
                        indices[j] -= 1
            else:
                poly = poly_res * MyPoly.INV[coeff]
                bullets[i]['base'] = poly
                bullets[i]['mon'] = mon
                if bullets[i]['type'] == TYPE_DIFF:  #
                    bullets[i]['base_diff'] = (bullets[i]['base_diff'] -
                                               sum(bullets[j]['base_diff'] * proj[j-index_nb]
                                                   for j in range(indices[IINDEX_DIFF], i))) * MyPoly.INV[coeff]
                elif bullets[i]['type'] == TYPE_IMAGE:
                    bullets[i]['base_inv_diff'] = (bullets[i]['base_inv_diff'] -
                                                   sum(bullets[j]['base_inv_diff'] * proj[j-index_nb]
                                                       for j in range(indices[IINDEX_IMAGE], i))) * MyPoly.INV[coeff]
                i += 1
        return result

    def _gen_basis_new_page(self):
        """ executed when a new bullet is added """
        for deg, bullets in self.entries[-1].items():
            indices = self.indices[-1][deg]
            for i in range(indices[IINDEX_IMAGE], len(bullets)):
                poly_res = poly_mod(bullets[i]['poly'], bullets[:i])
                mon, coeff = poly_res.get_item()
                poly = poly_res * MyPoly.INV[coeff]
                bullets[i]['base'] = poly
                bullets[i]['mon'] = mon

    def add_bullet(self, bullet: dict, deg: Optional[tuple] = None) -> int:
        if deg is None:
            deg = bullet['poly'].deg()
        if self.page > 2 and bullet['type'] == TYPE_REL:
            raise SSValueError("{}: Please add relations only on page 2".format(deg))
        if deg not in self.entries[-1]:
            if bullet['type'] == TYPE_TBD:
                self.indices[-1][deg] = [0, 0, 0]
            elif bullet['type'] == TYPE_REL:
                self.indices[-1][deg] = [1, 1, 1]
            else:
                raise SSValueError("{}: Illegal new bullet".format(deg))
            self.entries[-1][deg] = [bullet]
            return NB_NEW
        else:
            iindex = IINDEX_INSERT[bullet['type']]
            if iindex is not None:
                insert_index = self.indices[-1][deg][iindex]
                self.entries[-1][deg].insert(insert_index, bullet)
                for i in range(iindex, IINDEX_NUM):
                    self.indices[-1][deg][i] += 1
                if bullet['base'] is None:
                    result = self._gen_basis_new_bullet(deg, insert_index)
                    return result
                return NB_NEW
            else:
                self.entries[-1][deg].append(bullet)
                return NB_NEW

    def add_single_gen(self, g: str, deg: tuple):
        """ add a single bullet """
        if self.multiplicative:
            raise SSClassError("This spectral sequence requires multiplicative structure")
        if not self._inside_range(deg):
            print("generator outside of the range")
            return
        gen_poly = MyPoly.gen(g, deg)
        bullet = {'poly': gen_poly, 'base': gen_poly,
                  'mon': ((g, 1),), 'diff': None, 'type': TYPE_TBD}
        self.add_bullet(bullet, deg)

    def add_free_gen(self, g: str, deg: tuple) -> Optional[MyPoly]:
        """ add a generator g with deg freely """
        if not self.multiplicative:
            raise SSClassError("This spectral sequence does not support multiplicative structure")
        if not self._inside_range(deg):
            print("generator outside of the range")
            return None
        gen_poly = MyPoly.gen(g, deg)
        for i, j in self._get_degs("bullet", deg):
            for bullet in self.entries[0][(i, j)]:
                bul = {'poly': gen_poly * bullet['poly'], 'base': gen_poly * bullet['base'],
                       'mon': MyPoly.prod_mons(((g, 1),), bullet['mon']),
                       'diff': None, 'type': TYPE_REL if bullet['type'] == TYPE_REL else TYPE_TBD}
                self.add_bullet(bul, (deg[0] + i, deg[1] + j))
        return gen_poly

    def add_diff(self, poly_saddr1, poly_saddr2, b_log=True):
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
        self.x_max.append(self.x_max[-1] - self.page)
        self.y_max.append(self.y_max[-1] - self.page + 1)
        entries_page = dict((deg, self.homology(deg)) for deg in self.entries[-1]
                            if deg[0] <= self.x_max[-1] and deg[1] <= self.y_max[-1])
        indices_page = dict((deg, self.gen_indices(bullets)) for deg, bullets in entries_page.items())
        self.entries.append(entries_page)
        self.indices.append(indices_page)
        self.page += 1
        self._gen_basis_new_page()
        self.deduce_init()
