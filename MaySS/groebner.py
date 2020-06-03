"""Algebras based on Groebner basis.

Tailored for May spectral sequence.
The algebra is tri-graded with (s, t, u).
Monomials are structured by sparse vectors.
The Groebner basis is truncated by a predicate function."""
# TODO: use sqlite to store basis
# TODO: consider square-free property of E_2
# TODO: special treatment of relations x_ix_j=0
# TODO: support json
import copy
import heapq
import pickle
from bisect import bisect_left
from collections import defaultdict
from itertools import chain, repeat, combinations, groupby, combinations_with_replacement
from typing import Tuple, List, Dict, Set, Type, Iterable, NamedTuple, Union, Callable, Any
from algebras import BaseAlgebras as BA, linalg
from algebras.mymath import add_dtuple, sub_dtuple, div_mod_dtuple, Vector, get_one_element, \
    tex_parenthesis, tex_pow, two_expansion
_S_MAX = 1000


class Gen(NamedTuple):
    index: int
    name: str
    deg: int
    deg3d: Vector


class DgaGen(NamedTuple):
    index: int
    name: str
    deg: int
    deg3d: Vector
    diff: set


class RelCache(NamedTuple):
    deg: int
    is_rel_gen: bool
    rel: set
    deg3d: Vector = Vector((0, 0, 0))


class GbAlgMod2(BA.AlgebraMod2):
    """A factory for algebras using Groebner basis.

    `GbAlgMod2.new_alg()` creates a new algebra with
    its own generators and relations.
    """

    # generators: list of (index, name, degree)
    generators = None  # type: List[Gen]
    rels = None  # type: Dict[tuple, Set[Tuple[int, int]]]
    _rels_gen_leads = None  # type: Set[tuple]
    _rels_cache = None  # type: List[RelCache]
    key = None  # type: Callable[[tuple], Any]
    pred = None  # type: Callable[[tuple], Any]
    auto_simplify = None  # type: bool
    _attributes = ["generators", "rels", "_rels_gen_l"
                                         "eads", "_rels_cache", "key", "pred", "auto_simplify"]
    _name_index = 0

    @staticmethod
    def new_alg(*, key=None, pred=None) -> "Type[GbAlgMod2]":
        """Return a dynamically created subclass of GbAlgMod2.

        When key=None, use reversed lexicographical ordering by default."""
        cls = GbAlgMod2
        class_name = f"GbAlgMod2_{cls._name_index}"
        cls._name_index += 1
        dct = {'generators': [], 'rels': {}, '_rels_gen_leads': set(), '_rels_cache': [],
               'key': key, 'pred': pred or cls.pred_default, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbAlgMod2]":
        """Return a copy of current algebra."""
        class_name = f"GbAlgMod2_{GbAlgMod2._name_index}"
        GbAlgMod2._name_index += 1
        dct = {'generators': cls.generators.copy(), 'rels': copy.deepcopy(cls.rels),
               '_rels_gen_leads': cls._rels_gen_leads.copy(),
               '_rels_cache': copy.deepcopy(cls._rels_cache),
               'key': cls.key, 'pred': cls.pred, 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbAlgMod2,), dct)

    @classmethod
    def save_alg(cls, filename):
        """Save to a pickle file."""
        with open(filename, 'wb') as file:
            pickle.dump(["5-11-2020"] + [getattr(cls, attr) for attr in cls._attributes], file)

    @classmethod
    def load_alg(cls, filename) -> Union[Type["GbAlgMod2"], Type["GbDga"]]:
        """Create an algebra from a pickle file."""
        with open(filename, 'rb') as file:
            init_list = pickle.load(file)
            if cls is GbAlgMod2:
                class_name = f"GbAlgMod2_{cls._name_index}"
            elif cls is GbDga:
                class_name = f"GbDga_{cls._name_index}"
            else:
                raise BA.MyKeyError("cls should be either GbAlgMod2 or GbDga.")
            cls._name_index += 1
            if init_list[0] == "5-11-2020":
                dct = {attr: init_v for attr, init_v in zip(cls._attributes, init_list[1:])}
            else:
                raise ValueError("file version not recognized")
            # noinspection PyTypeChecker
            return type(class_name, (cls,), dct)

    @staticmethod
    def construct_alg(pred, basis, cls_basis, fn: Callable[[Any], set], naming: Callable):
        """Construct an algebra from basis."""
        R = GbAlgMod2.new_alg(pred=pred)
        # noinspection PyArgumentList
        R_basis_mons = defaultdict(list, {Vector((0, 0, 0)): [()]})
        map_alg = linalg.GradedLinearMapKMod2()
        image_gens = {}
        basis_V = basis
        if basis:
            d = get_one_element(basis)
            if type(basis[d]) is not linalg.VectorSpaceMod2:
                basis_V = {d: linalg.VectorSpaceMod2({m} for m in basis[d]) for d in basis}

        ds = sorted(d for d in basis if pred(d))
        for d in ds:
            if d == (0, 0, 0):
                continue
            BA.Monitor.print(f"{d=}")
            for x in (basis_V[d] / map_alg.image(d)).basis(cls_basis):
                BA.Monitor.print(f"{x=}", 1)
                gen_name = naming(x, d)
                R.add_gen(gen_name, d[1], d)
                index = R.generators[-1][0]
                image_gens[gen_name] = x

                ds_R_basis_mons = {d: len(R_basis_mons[d]) for d in R_basis_mons}
                leadings = []
                for d1 in ds_R_basis_mons:
                    if pred(d1 + d):
                        for i in range(ds_R_basis_mons[d1]):
                            m1 = R_basis_mons[d1][i]
                            e = 1
                            while pred(d2 := d1 + d * e):
                                m2 = m1 + ((index, -e),)  # type: tuple
                                if any(map(le_dtuple, leadings, repeat(m2))):
                                    break
                                R_basis_mons[d2].append(m2)
                                r2 = R(m2)
                                fr2 = fn(r2.evaluation(image_gens))
                                map_alg.add_maps_set([(r2.data, fr2)], d2)
                                if map_alg.kernel(d2):
                                    R.add_rels_data(map_alg.kernel(d2).basis(set))
                                    map_alg.kernel(d2).clear()
                                    leadings = [m for m in R.rels if m[-1][0] == index]
                                    break
                                e += 1
        R.add_rels_cache()
        return R, image_gens, map_alg

    # ----- AlgebraMod2
    @classmethod
    def mul_mons(cls, mon1: tuple, mon2: tuple):
        m = add_dtuple(mon1, mon2)
        return cls.simplify_data({m}) if cls.auto_simplify else m

    @classmethod
    def str_mon(cls, mon: tuple):
        if mon:
            return "".join(tex_pow(cls.get_gen(i).name, -e) for i, e in mon)
        else:
            return "1"

    @classmethod
    def repr_mon(cls, mon: tuple, clsname):
        if mon:
            return " * ".join(f"{clsname}.gen(\"{cls.get_gen(i).name}\") ** {-e}"
                              if -e > 1 else f"{clsname}.gen(\"{cls.get_gen(i).name}\")"
                              for i, e in mon)
        else:
            return f"{clsname}.unit()"

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(cls.get_gen(i).deg * -e for i, e in mon)

    @classmethod
    def deg3d_mon(cls, mon: tuple):
        return sum((cls.get_gen(i).deg3d * -e for i, e in mon), Vector((0, 0, 0)))

    def deg(self):
        """Require `self` to be homogeneous."""
        for m in self.data:
            return self.deg_mon(m)

    def deg3d(self):
        """Require `self` to be homogeneous."""
        for m in self.data:
            return self.deg3d_mon(m)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, name: str, deg, deg3d=(0, 0, 0)):
        """Add a new generator and return it."""
        if cls.pred(deg3d):
            index = cls.generators[-1][0] + 1 if cls.generators else 0
            cls.generators.append(Gen(index, name, deg, Vector(deg3d)))
            m = ((index, -1),)
            return cls(m).simplify() if cls.auto_simplify else cls(m)

    @classmethod
    def add_gens(cls, name_deg_deg3d_s):
        """Add generators. name_deg_deg3d_s is a list of tuples (name, deg, deg3d)."""
        index = cls.generators[-1][0] + 1 if cls.generators else 0
        for n, d, d3d in name_deg_deg3d_s:
            if cls.pred(d3d):
                cls.generators.append(Gen(index, n, d, Vector(d3d)))
                index += 1

    @classmethod
    def remove_gen(cls, name: str):
        """Remove a generator.

        If the generator `name` equals zero in the algebra whose relations are simplified,
        call this function to remove this generator. Use only when _rels_cache is empty."""
        for i, item in enumerate(cls.generators):
            if item[1] == name:
                break
        else:
            raise ValueError(f"no generator named {name}")
        m = ((item[0], -1),)
        assert not cls.rels[m]
        del cls.rels[m]
        del cls.generators[i]

    @classmethod
    def rename_gen(cls, old_name, new_name):
        """Rename a generator."""
        for i, gen in enumerate(cls.generators):
            if gen.name == old_name:
                break
        else:
            raise ValueError(f"no generator named {old_name}")
        cls.generators[i] = Gen(gen.index, new_name, gen.deg, gen.deg3d)

    @classmethod
    def reorder_gens(cls, index_map=None, key=None):  # TODO: create a new alg instead
        """Reorganize the relations by a new ordering of generators and a new key function.
        The new i'th generator is the old `index_map[i]`'th generator."""
        index_map_inv = {}
        for i, fi in enumerate(index_map):
            index_map_inv[fi] = i
        num_gens = len(cls.generators)
        rel_generators = cls.get_rel_gens()
        cls.key = key
        if index_map:
            def f(m):
                return tuple((index_map_inv[_i], _e) for _i, _e in m)
            assert num_gens == len(index_map)
            rel_generators = [{f(m) for m in rel} for rel in rel_generators]
            cls.generators = [Gen(i, (gen := cls.get_gen(index_map[i])).name, gen.deg, gen.deg3d)
                              for i in range(num_gens)]
        cls.rels = {}
        cls._rels_gen_leads = set()
        cls.add_rels_data(rel_generators, clear_cache=True)

    @classmethod
    def add_rel_data(cls, rel: set, clear_cache=False):
        """Add relations."""
        if rel:
            deg3d = cls.deg3d_data(rel)
            if cls.pred(deg3d):
                deg = cls.deg_data(rel)
                heapq.heappush(cls._rels_cache, RelCache(deg, True, rel, deg3d))
                cls.add_rels_cache(None if clear_cache else deg)

    @classmethod
    def add_rels_data(cls, rels: Iterable[set], sorted_=False, clear_cache=False):
        """Add relations."""
        for rel in (rels if sorted_ else sorted(rels, key=cls.deg_data)):
            cls.add_rel_data(rel)
        if clear_cache:
            cls.add_rels_cache()

    @classmethod
    def add_rel(cls, rel: "GbAlgMod2", clear_cache=False):
        """Add a relation."""
        if not rel.is_homo():
            raise ValueError(f'relation {rel} not homogeneous!')
        cls.add_rel_data(rel.data, clear_cache)

    @classmethod
    def add_rels(cls, rels: Iterable["GbAlgMod2"], sorted_=False, clear_cache=False):
        """Add a relation."""
        for rel in rels:
            if not rel.is_homo():
                raise ValueError(f'relation {rel} not homogeneous!')
        cls.add_rels_data((rel.data for rel in rels), sorted_, clear_cache)

    @classmethod
    def simplify_rels(cls):
        """Simplify `cls.rels`."""
        for m in cls.rels:
            cls.rels[m] = cls.simplify_data(cls.rels[m])
        for i, rel_cache in enumerate(cls._rels_cache):
            cls._rels_cache[i] = RelCache(rel_cache.deg, rel_cache.is_rel_gen, cls.simplify_data(rel_cache.rel))

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    # getters --------------------------
    @classmethod
    def get_gen(cls, index: int) -> Union[Gen, DgaGen]:
        i = bisect_left(cls.generators, (index, "", 0))
        return cls.generators[i]

    @classmethod
    def get_gen_by_name(cls, name: str) -> Union[Gen, DgaGen]:
        for gen in cls.generators:
            if gen.name == name:
                return gen
        raise BA.MyKeyError(f"generator {name} not found")

    @classmethod
    def deg_data(cls, data: set):
        """Return the degree of `data`."""
        for m in data:
            return cls.deg_mon(m)

    @classmethod
    def deg3d_data(cls, data: set):
        """Return the degree of `data`."""
        for m in data:
            return cls.deg3d_mon(m)

    @classmethod
    def get_lead(cls, data):
        """Return the leading term of `data`."""
        return max(data, key=cls.key) if cls.key else max(data)

    @classmethod
    def get_num_gens(cls):
        """Return the number of generators."""
        return len(cls.generators)

    @classmethod
    def get_cache_size(cls):
        return len(cls._rels_cache)

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        for gen in cls.generators:
            if gen.name == k:
                m = ((gen.index, -1),)
                return cls(m).simplify() if cls.auto_simplify else cls(m)
        else:
            raise BA.MyKeyError(f"No generator named {k}")

    @classmethod
    def get_rel_gens(cls):
        """Return a minimal generating set of `cls.rels` or `ideal`, ordered by degree."""
        rel_gens = []
        for m in sorted(cls._rels_gen_leads, key=cls.deg_mon):
            rel = cls.rels[m] | {m}
            rel_gens.append(rel)
        return rel_gens

    @classmethod
    def get_ideal_gens(cls, ideal: List[set]):
        """Return the minimal generating set of `ideal`."""
        A = cls.copy_alg()
        rel_gens = []
        for rel in sorted(ideal, key=A.deg_data):
            A.add_rels_cache(A.deg_data(rel))
            rel = A.simplify_data(rel)
            if rel:
                rel_gens.append(rel)
                A.add_rel_data(rel)
        return rel_gens

    @classmethod
    def basis_mons(cls, pred=None, basis: Dict[Vector, List[tuple]] = None):
        """Return a list of basis (mon, deg)."""
        pred = pred or cls.pred
        # noinspection PyArgumentList
        result = basis or defaultdict(list, {Vector((0, 0, 0)): [()]})
        old_ds = set(result)
        leadings = sorted(cls.rels, key=lambda _m: _m[-1][0])
        leadings = {index: list(g) for index, g in groupby(leadings, key=lambda _m: _m[-1][0])}
        for gen in cls.generators:
            index, deg3d = gen.index, gen.deg3d
            ds = list(result)
            for d in ds:
                if (d_ := d + deg3d) not in old_ds and pred(d_):
                    for m in result[d]:
                        i_m = m[-1][0] if m else -1
                        if index == i_m and d in old_ds:
                            e = 1
                            while pred(d1 := d + deg3d * e):
                                m1 = m[:-1] + ((m[-1][0], m[-1][1] - e),)
                                if index in leadings and any(map(le_dtuple, leadings[index], repeat(m1))):

                                    break
                                result[d1].append(m1)
                                e += 1
                        elif index > i_m:
                            e = 1
                            while pred(d1 := d + deg3d * e):
                                m1 = m + ((index, -e),)
                                if index in leadings and any(map(le_dtuple, leadings[index], repeat(m1))):
                                    break
                                elif d1 in result:
                                    result[d1].append(m1)
                                else:
                                    result[d1] = [m1]
                                e += 1
        return result

    @classmethod
    def basis(cls, pred=None):
        if type(pred) is int:
            basis = cls.basis_mons(pred=lambda d3d: d3d[0] <= pred)
        else:
            basis = cls.basis_mons(pred)
        for basis_d in basis.values():
            for m in basis_d:
                yield cls(m)

    @classmethod
    def basis_mons_h0(cls, pred=None):
        """Return a list of $h_0$-structure lines."""
        pred = pred or cls.pred

        def d3d_h0(d3d):
            return d3d - ((d3d0 := d3d[0]), d3d0, d3d0)

        # noinspection PyArgumentList
        result = defaultdict(list, {Vector((0, 0, 0)): [((0, -_S_MAX),)]})
        old_ds = set(result)
        leadings = sorted(cls.rels, key=lambda _m: _m[-1][0])
        leadings = {index: list(g) for index, g in groupby(leadings, key=lambda _m: _m[-1][0])}
        for gen in cls.generators:
            if gen.index == 0:
                continue
            index, deg3d = gen.index, d3d_h0(gen.deg3d)
            ds = list(result)
            for d in ds:
                if (d_ := d + deg3d) not in old_ds and pred(d_):
                    for m in result[d]:
                        i_m = m[-1][0] if m else -1
                        if index == i_m and d in old_ds:
                            e = 1
                            while pred(d1 := d + deg3d * e):
                                m1 = m[:-1] + ((m[-1][0], m[-1][1] - e),)
                                if index in leadings and any(map(le_dtuple, leadings[index], repeat(m1))):
                                    break
                                result[d1].append(m1)
                                e += 1
                        elif index > i_m:
                            e = 1
                            while pred(d1 := d + deg3d * e):
                                m1 = m + ((index, -e),)
                                if index in leadings and any(map(le_dtuple, leadings[index], repeat(m1))):
                                    break
                                result[d1].append(m1)
                                e += 1
        return result

    @classmethod
    def invert_h0(cls, pred=None):
        """Return the localization with h_0 inverted."""
        pred = pred or cls.pred
        basis = cls.basis_mons_h0(pred)

        def fn(x):
            return {((0, -_S_MAX),) + _m[1:] for _m in x.data}

        def naming(x, _):
            y = type(x)({_m[1:] for _m in x.data})
            return str(y)
        return cls.construct_alg(pred, basis, cls, fn, naming)[0]

    @classmethod
    def is_reducible(cls, mon):
        """Determine if mon is reducible by `cls.rels`."""
        return any(le_dtuple(m, mon) for m in cls.rels)

    def evaluation(self, image_gens: Dict[str, "GbAlgMod2"]):
        """Return f(self) where f is an algebraic map determined by `image_gens`."""
        for v in image_gens.values():
            R = type(v)
            break
        else:
            raise ValueError("empty image_gens")
        zero = R.zero() if issubclass(R, BA.Algebra) else 0
        unit = R.unit() if issubclass(R, BA.Algebra) else 1
        result = zero
        for m in self.data:
            fm = unit
            for i, e in m:
                fm *= image_gens[self.get_gen(i).name] ** (-e)
            result += fm
        return result

    @classmethod
    def print_latex_alg(cls, show_gb=True):
        """For latex."""
        print("\\section{Generators}\n")
        # noinspection PyArgumentList
        gens = defaultdict(list)
        for item in cls.generators:
            gens[item.deg3d[0]].append(item)
        for s in sorted(gens):
            print(f"\\subsection*{{s={s}}}\n")
            print(', '.join(f'${item.name}$' for item in gens[s]), end="\\vspace{3pt}\n\n")
        print("\\section{Relations}\n")
        if show_gb:
            for m in cls.rels:
                if m in cls._rels_gen_leads:
                    print(f"${cls(m)} = {cls(cls.rels[m])}$\\vspace{{3pt}}\n")
                else:
                    print(f"$\\bullet\\hspace{{4pt}} {cls(m)} = {cls(cls.rels[m])}$\\vspace{{3pt}}\n")
        else:
            for m in cls._rels_gen_leads:
                print(f"${cls(m)} = {cls(cls.rels[m])}$\n")

    @classmethod
    def markdown_alg(cls, show_all=True):
        """For Jupyter notebook."""
        from IPython.display import Markdown
        result = f"#### Generators:\n${', '.join(item[1] for item in cls.generators)}$\n\n"
        result += "#### Relations:\n"

        if show_all:
            result += "\\begin{align*}"
            for m in sorted(cls.rels, key=cls.deg_mon):
                if m in cls._rels_gen_leads:
                    result += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
                else:
                    result += f"\\bullet\\hspace{{4pt}} {cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            result += "\\end{align*}\n\n"
        else:
            result += "\\begin{align*}"
            for m in sorted(cls._rels_gen_leads, key=cls.deg_mon):
                result += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            result += "\\end{align*}\n\n"
        return Markdown(result)

    @classmethod
    def latex_ideal(cls, gb: Iterable[set]):
        from IPython.display import Latex
        result = "\\begin{align*}"
        for data in gb:
            result += f"&{cls(data)}\\\\\n"
        result += "\\end{align*}"
        return Latex(result)

    @staticmethod
    def latex_annilators(annilators:  List[List[Tuple["GbAlgMod2", str, int]]]):
        from IPython.display import Latex
        result = "\\begin{align*}\n"
        for a in annilators:
            s = "+".join(f"{name}{tex_parenthesis(coeff)}" for coeff, name, d in a)
            result += f"& {s}=0\\\\\n"
        result += "\\end{align*}"
        return Latex(result)

    @staticmethod
    def repr_annilators(clsname, annilators:  List[List[Tuple["GbAlgMod2", str, int]]]):
        result = ""
        for a in annilators:
            s = " + ".join(f"{clsname}.gen(\"{name}\") * {tex_parenthesis(coeff.repr_(clsname))}"
                           for coeff, name, d in a)
            result += f"{s}\n"
        return result

    @classmethod
    def to_dga(cls, deg_diff: tuple) -> "Type[GbDga]":
        class_name = f"GbDGA_{GbDga._name_index}"
        GbDga._name_index += 1
        # noinspection PyTypeChecker
        dct = {'generators': [DgaGen(*gen, None) for gen in cls.generators], 'rels': copy.deepcopy(cls.rels),
               '_rels_gen_leads': cls._rels_gen_leads.copy(),
               '_rels_cache': copy.deepcopy(cls._rels_cache),
               'key': cls.key, 'pred': cls.pred, 'auto_simplify': cls.auto_simplify,
               'deg_diff': Vector(deg_diff)}
        # noinspection PyTypeChecker
        return type(class_name, (GbDga,), dct)

    @classmethod
    def pred_default(cls, _) -> bool:
        return True

    def is_gen(self):
        if len(self.data) == 1:
            for m in self.data:
                if sum(-e for i, e in m) == 1:
                    return True
        return False

    # algorithms --------------------------
    @classmethod
    def simplify_data(cls, data: set):
        """Return simplified `data`. `data` will remain unchanged."""
        state_backup = cls.auto_simplify
        cls.auto_simplify = False
        s = data.copy()
        result = set()
        while s:
            mon = cls.get_lead(s)
            s.remove(mon)
            for m in cls.rels:
                if le_dtuple(m, mon):
                    q, r = div_mod_dtuple(mon, m)
                    m_to_q = (cls(cls.rels[m]) ** q).data
                    s ^= {add_dtuple(r, m1) for m1 in m_to_q}
                    break
            else:
                result ^= {mon}
        cls.auto_simplify = state_backup
        return result

    @classmethod
    def add_rels_cache(cls, deg_max=None):
        """Add relations from cache up to deg `deg_max`."""
        hq = cls._rels_cache
        while hq and (deg_max is None or hq[0][0] <= deg_max):
            deg, is_rel_gen, r, _ = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                m = cls.get_lead(r)
                if is_rel_gen:
                    cls._rels_gen_leads.add(m)
                redundant_leading_terms = []
                for m1, v1 in cls.rels.items():
                    if gcd_nonzero_dtuple(m, m1):
                        lcm = max_dtuple(m, m1)
                        if cls.pred(cls.deg3d_mon(lcm)):
                            dif = sub_dtuple(lcm, m)
                            dif1 = sub_dtuple(lcm, m1)
                            new_rel = {add_dtuple(_m, dif) for _m in r}
                            v1dif1 = {add_dtuple(_m, dif1) for _m in v1}
                            new_rel -= {lcm}
                            new_rel ^= v1dif1
                            if le_dtuple(m, m1):
                                BA.Monitor.count("redundant_leading_terms")  #
                                redundant_leading_terms.append(m1)
                                is_lead = m1 in cls._rels_gen_leads
                                if new_rel:
                                    heapq.heappush(hq, RelCache(cls.deg_mon(lcm), is_lead, new_rel))
                                if is_lead:
                                    cls._rels_gen_leads.remove(m1)
                            elif new_rel:
                                heapq.heappush(hq, RelCache(cls.deg_mon(lcm), False, new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls.rels[m_redundant]
                cls.rels[m] = r - {m}

    @classmethod
    def get_vector_gens(cls, ideal: List[List[Tuple["GbAlgMod2", str, int]]], *, inplace=False):
        """Return the minimal generating set of `ideal`, which is a A-submodule of A^n.

        `ideal` is a list of [(ele, name, deg), ...].
        The names should not overlap with existing generator names of `cls`."""
        # TODO: create class AugModuleMod2
        A = cls if inplace else cls.copy_alg()
        rels = []
        added_names = set()
        for index, v in enumerate(ideal):
            rel = A.zero()
            for ele, name, deg in v:
                name = f"v_{{{name}}}"
                x = A.gen(name) if name in added_names else A.add_gen(name, deg)  #
                added_names.add(name)
                rel += x * ele
            if rel:
                rels.append((index, rel))
        rels_module = []
        for i1, i2 in combinations_with_replacement(added_names, 2):
            rels_module.append(A.gen(i1) * A.gen(i2))
        A.add_rels(rels_module)
        result = []
        for index, rel in sorted(rels, key=lambda _x: _x[1].deg()):
            A.add_rels_cache(rel.deg())
            rel.simplify()
            if rel:
                result.append(index)
                A.add_rel(rel)
        return [ideal[i] for i in result]

    @classmethod
    def ann(cls, x):
        """Return the groebner basis for the ideal {y | xy=0}."""
        annihilators = cls.ann_seq([(x, '')])
        return [a[0][0].data for a in annihilators]

    @classmethod
    def ann_seq(cls, ele_names: List[Tuple["GbAlgMod2", str]]):
        """Return relations among elements: $\\sum a_ie_i=0$."""
        A = cls.copy_alg()
        index = A.generators[-1][0] if A.generators else -1
        if cls.key:
            A.key = lambda _m: (_m[bisect_left(_m, (index, 0)):], cls.key(_m))
        else:
            A.key = lambda _m: (_m[bisect_left(_m, (index, 0)):], _m)
        rels_new = []
        for ele, name in ele_names:
            x = A.add_gen(name, ele.deg())
            rels_new.append(x + ele)
        A.add_rels(rels_new, clear_cache=True)
        annilators = []
        for m in A.rels:
            if m[-1][0] > index:
                a = []
                for m1 in chain((m,), A.rels[m]):
                    gen = A.get_gen(m[-1][0])
                    m11 = m1[:-1] + (m1[-1][0], m1[-1][1] + 1) if m1[-1][1] + 1 != 0 else m1[:-1]
                    a.append((m11, gen.name, gen.deg))
                annilators.append(a)
        for en1, en2 in combinations(ele_names, 2):
            ele1, name1 = en1
            deg1 = ele1.deg()
            ele2, name2 = en2
            deg2 = ele2.deg()
            a = []
            for m1 in ele1.data:
                a.append((m1, name2, deg2))
            for m2 in ele2.data:
                a.append((m2, name1, deg1))
            annilators.append(a)
        if cls.key:
            def key(_m):
                return A.deg_mon(_m[bisect_left(_m, (index, 0)):]), cls.key(_m)
        else:
            def key(_m):
                return A.deg_mon(_m[bisect_left(_m, (index, 0)):]), _m
        A.reorder_gens(key=key)
        annilators = [[(cls(A.simplify_data({_m})), name, deg) for _m, name, deg in a] for a in annilators]
        return cls.get_vector_gens(annilators)

    @classmethod
    def subalgebra(cls, ele_names: List[Tuple["GbAlgMod2", str]], *, key=None):
        """Return the subalgebra generated by `ele_names`."""
        A = cls.copy_alg()
        index = A.generators[-1][0] if A.generators else -1

        def key1(_m):
            _m1, _m2 = _m[:(i := bisect_left(_m, (index, 0)))], _m[i:]
            return cls.deg_mon(_m1), (cls.key(_m1) if cls.key else _m1), (key(_m2) if key else _m2)
        A.key = key1
        for ele, name in ele_names:
            x = A.add_gen(name, ele.deg())
            A.add_rel(x + ele)
        A.add_rels_cache()
        A.generators = A.generators[bisect_left(A.generators, (index + 1,)):]
        A.key = key
        A._rels_gen_leads = set()
        rels, A.rels = A.rels, {}
        for m in rels:
            if m[0][0] > index:
                rel_subalg = {m[bisect_left(m, (index, 0)):]} | {_m[bisect_left(_m, (index, 0)):] for _m in rels[m]}
                A.add_rel_data(rel_subalg)
        A.add_rels_cache()
        return A


class GbDga(GbAlgMod2):
    """A factory for DGA over F_2."""
    generators = None  # type: List[DgaGen]
    deg_diff = None
    _attributes = ["generators", "rels", "_rels_gen_leads", "_rels_cache", "key", "pred", "auto_simplify", "deg_diff"]

    @staticmethod
    def new_alg(*, key=None, pred=None, deg_diff=None) -> "Type[GbDga]":
        """Return a dynamically created subclass of GbDga."""
        cls = GbAlgMod2
        class_name = f"GbAlgMod2_{cls._name_index}"
        cls._name_index += 1
        if deg_diff is not None:
            deg_diff = Vector(deg_diff)
        else:
            raise BA.MyDegreeError("degree of differential not supplied")
        dct = {'generators': [], 'rels': {}, '_rels_gen_leads': set(), '_rels_cache': [],
               'key': key, 'pred': pred, 'auto_simplify': True, 'deg_diff': deg_diff}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbDga]":
        """Return a copy of current algebra."""
        class_name = f"GbAlgMod2_{GbAlgMod2._name_index}"
        GbAlgMod2._name_index += 1
        dct = {'generators': cls.generators.copy(), 'rels': copy.deepcopy(cls.rels),
               '_rels_gen_leads': cls._rels_gen_leads.copy(),
               '_rels_cache': copy.deepcopy(cls._rels_cache),
               'key': cls.key, 'pred': cls.pred, 'auto_simplify': cls.auto_simplify,
               'deg_diff': cls.deg_diff}
        # noinspection PyTypeChecker
        return type(class_name, (GbAlgMod2,), dct)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, name: str, deg, deg3d=(0, 0, 0), diff=None):
        """Add a new generator and return it."""
        index = cls.generators[-1][0] + 1 if cls.generators else 0
        if diff is None:
            diff = set()
        elif type(diff) is not set:
            diff = diff.data
        if diff and cls.deg3d_data(diff) - deg3d != cls.deg_diff:
            raise BA.MyDegreeError("inconsistent differential degree")
        cls.generators.append(DgaGen(index, name, deg, Vector(deg3d), diff))
        m = ((index, -1),)
        return cls(m).simplify() if cls.auto_simplify else cls(m)

    @classmethod
    def set_diff(cls, gen_name: str, diff: Union[None, set, "GbDga"]):
        """Define the differential of gen_name."""
        for i, gen in enumerate(cls.generators):
            if gen.name == gen_name:
                break
        else:
            raise BA.MyKeyError(f"generator {gen_name} not found")
        if i is not None:
            if type(diff) is not set and diff is not None:
                diff = diff.data
            if diff and (not cls(diff).is_homo() or cls.deg3d_data(diff) - gen.deg3d != cls.deg_diff):
                raise BA.MyDegreeError("inconsistent differential degree")
            gen = cls.generators[i]
            cls.generators[i] = DgaGen(gen.index, gen.name, gen.deg, gen.deg3d, diff)

    def diff(self):
        """Return the boundary of the chain."""
        result = set()
        for m in self.data:
            for (i, (index, e)) in enumerate(m):
                if e % 2:
                    m1 = m[:i] + ((index, e + 1),) + m[i+1:] if e + 1 else m[:i] + m[i+1:]
                    m1_by_dg_i = {add_dtuple(m1, _m) for _m in self.get_gen(index).diff}
                    result ^= m1_by_dg_i
        return type(self)(result).simplify()

    @classmethod
    def rename_gen(cls, old_name, new_name):
        """Rename a generator."""
        for i, gen in enumerate(cls.generators):
            if gen.name == old_name:
                break
        else:
            raise ValueError(f"no generator named {old_name}")
        cls.generators[i] = DgaGen(gen.index, new_name, gen.deg, gen.deg3d, gen.diff)

    # getters ----------------------------
    @classmethod
    def homology(cls, pred, BH=None, page=4) -> Tuple[Type[GbAlgMod2], dict, linalg.GradedLinearMapKMod2]:
        """Compute HA. Return (HA, list of representing cycles)."""
        B, H = BH or cls.basis_BH(pred)
        gen_index = [0]

        def fn(x):
            return B[d3d].res_set(x.data) if (d3d := x.deg3d()) in B else set()

        def naming(x, d):
            if x.is_gen():
                gen_name = str(x)
            elif (name := str(x)).count("+") <= 1:
                gen_name = f'[{name}]'
            else:
                gen_name = f"x_{{{page}, {d[1] - d[0]}, {gen_index[0]}}}"
                gen_index[0] += 1
            return gen_name
        R, image_gens, map_alg = cls.construct_alg(pred, H, cls, fn, naming)
        return R, image_gens, map_alg

    @classmethod
    def homology_h0(cls, pred, BH, page) -> Tuple[Type[GbAlgMod2], dict, linalg.GradedLinearMapKMod2]:
        """Compute h_0^{-1}HA. Return (HA, list of representing cycles, map for localization)."""
        B, H = BH or cls.basis_BH_h0(pred)
        gen_index = defaultdict(int)

        def fn_invert_h0(x):
            data = set()
            for m in x.data:
                data ^= {((0, -_S_MAX),) + m[1:]}
            d3d = cls.deg3d_data(data)
            if d3d:
                s = d3d[0]
                d3d = d3d - (s, s, s)
            return d3d, data

        def fn(x):
            d3d, data = fn_invert_h0(x)
            return B[d3d].res_set(data) if d3d in B else set()

        def naming(x, d):
            y = type(x)({_m[1:] for _m in x.data})
            if y.is_gen():
                gen_name = str(y)
            elif (name := str(y)).count("+") <= 1:
                gen_name = f'[{name}]'
            else:
                if page is None:
                    gen_name = f"x_{{{d[1]}, {gen_index[d[1]]}}}"
                else:
                    gen_name = f"x_{{{page}, {d[1]}, {gen_index[d[1]]}}}"
                gen_index[d[1]] += 1
            return gen_name
        R, image_gens, map_alg = cls.construct_alg(pred, H, cls, fn, naming)
        return R, image_gens, map_alg

    @classmethod
    def basis_BH(cls, pred, basis=None, BH=None):
        """Return the basis of cycles and homologies."""
        basis = basis or cls.basis_mons(pred=lambda _d: pred(_d) or pred(_d + cls.deg_diff))
        B, H = BH or ({}, {})
        map_diff = linalg.GradedLinearMapKMod2()
        print_end = " " * 10 + "\r"
        for d in sorted(basis):
            print(d, end=print_end)
            if pred(d) or pred(d + cls.deg_diff):
                if d not in B or d + cls.deg_diff not in B:
                    B_d = set(map_diff.image(d - cls.deg_diff).get_mons())
                    map_diff.add_maps_set((((r := cls(m)).data, r.diff().data) for m in basis[d]
                                           if m not in B_d), d)
        ds = sorted(d for d in basis if pred(d) and d not in B)
        B.update((d, map_diff.image(d - cls.deg_diff)) for d in ds)
        H.update((d, map_diff.kernel(d)) for d in ds)
        return B, H

    @classmethod
    def basis_BH_h0(cls, pred, basis=None, BH=None):
        """Return the basis of cycles and homologies with $h_0$ inverted."""
        s = cls.deg_diff[0]
        deg_diff = cls.deg_diff - (s, s, s)
        basis = basis or cls.basis_mons_h0(pred=lambda _d: pred(_d) or pred(_d + deg_diff))
        B, H = BH or ({}, {})
        map_diff = linalg.GradedLinearMapKMod2()
        print_end = " " * 10 + "\r"

        def fn_invert_h0(x):
            data = set()
            for m in x.data:
                data ^= {((0, -_S_MAX),) + m[1:]}
            return data
        for d in sorted(basis, key=lambda _d: (-_d[2], _d[1])):
            print(d, end=print_end)
            if pred(d) or pred(d + deg_diff):
                if d not in B or d + deg_diff not in B:
                    B_d = set(map_diff.image(d - deg_diff).get_mons())
                    map_diff.add_maps_set((((r := cls(m)).data, fn_invert_h0(r.diff())) for m in basis[d]
                                           if m not in B_d), d)
        ds = sorted(d for d in basis if pred(d) and d not in B)
        B.update((d, map_diff.image(d - deg_diff)) for d in ds)
        H.update((d, map_diff.kernel(d)) for d in ds)
        return B, H

    @classmethod
    def is_differentiable_mon(cls, mon):
        for i, e in mon:
            if cls.get_gen(i).diff is None:
                return False
        return True

    @classmethod
    def is_differentiable_data(cls, data):
        return all(map(cls.is_differentiable_mon, data))

    def is_differentiable(self):
        return self.is_differentiable_data(self.data)

    @staticmethod
    def contains_gen_mon(mon, index):
        """Return if the monomial contains the generator identified by index."""
        for i, e in mon:
            if i == index:
                return True
        return False

    @staticmethod
    def contains_gen_data(data, index):
        """Return if the data contains the generator identified by index."""
        return any(map(GbDga.contains_gen_mon, data, repeat(index)))

    @classmethod
    def determine_diff(cls, g: Union[str, int], basis: dict, image_gens=None):
        """Determine differentials by relations."""
        if type(g) is str:
            for gen in cls.generators:
                if gen.name == g:
                    break
            else:
                raise BA.MyKeyError(f"generator {g} not found")
        else:
            gen = cls.get_gen(g)
            g = gen.name
        deg_target = gen.deg3d + cls.deg_diff
        if deg_target not in basis:
            cls.set_diff(g, set())
            print(f"set d({g})=0")
            return
        # print("Possible summands:")
        # for m in basis[deg_target]:
        #     print(cls(m))
        rels = []
        cls.set_diff(g, set())
        for m in cls.rels:
            if cls.is_differentiable_mon(m) and all(map(cls.is_differentiable_mon, cls.rels[m])):
                if cls.contains_gen_mon(m, gen.index) or any(map(cls.contains_gen_mon, cls.rels[m], repeat(gen.index))):
                    rels.append({m} | cls.rels[m])
        possible_diffs = []
        for n in range(1 << len(basis[deg_target])):
            data = {basis[deg_target][i] for i in two_expansion(n)}
            if all(map(cls.is_differentiable_mon, data)) and cls(data).diff():
                continue
            cls.set_diff(g, data)
            compatible = True
            for rel in rels:
                if cls(rel).diff():
                    compatible = False
                    break
            if image_gens and image_gens[gen.name].is_differentiable() and\
                    cls(data).evaluation(image_gens) != image_gens[gen.name].diff():
                compatible = False
            if compatible:
                possible_diffs.append(data)
        if len(possible_diffs) == 1:
            cls.set_diff(g, possible_diffs[0])
            print(f"set d({g})={cls(possible_diffs[0])}")
        elif len(possible_diffs) == 0:
            raise BA.MyClassError(f"Invalid DGA. d({g})=?")
        else:
            for data in possible_diffs:
                print(f"d({g})={cls(data)} is possible.")
            cls.set_diff(g, None)

    @classmethod
    def determine_diffs(cls, basis: dict, image_gens=None):
        for gen in sorted((gen for gen in cls.generators if gen.diff is None),
                          key=lambda _gen: _gen.deg3d[1] - _gen.deg3d[0]):
            if gen.diff is None:
                cls.determine_diff(gen.index, basis, image_gens)

    @classmethod
    def print_latex_alg(cls, show_gb=False):
        """For latex."""
        super().print_latex_alg(show_gb)
        print("\\section{Differentials}\n")
        for gen in cls.generators:
            print(f"$d({gen.name})={cls(gen.diff)}$\\vspace{{3pt}}\n")

    @classmethod
    def markdown_alg(cls, show_all=False):
        """For Jupyter notebook."""
        from IPython.display import Markdown
        result = f"#### Generators:\n${', '.join(item[1] for item in cls.generators)}$\n\n"
        result += "#### Relations:\n"

        if show_all:
            result += "\\begin{align*}"
            for m in sorted(cls.rels, key=cls.deg_mon):
                if m in cls._rels_gen_leads:
                    result += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
                else:
                    result += f"\\bullet\\hspace{{4pt}} {cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            result += "\\end{align*}\n\n"
        else:
            result += "\\begin{align*}"
            for m in sorted(cls._rels_gen_leads, key=cls.deg_mon):
                result += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            result += "\\end{align*}\n\n"
        result += "#### Differentials:\n"
        result += "\\begin{align*}"
        for gen in cls.generators:
            result += f"d({gen.name}) &= {cls(gen.diff) if gen.diff is not None else '?'}\\\\\n"
        result += "\\end{align*}\n\n"
        return Markdown(result)


# operations for monomials with negative exponents
def le_dtuple(d1, d2):
    """Return if d1_i <= d2_i as sparse vectors."""
    d2_dict = dict(d2)
    return all(gen in d2_dict and exp >= d2_dict[gen] for gen, exp in d1)


def min_dtuple(d1, d2):
    """return (min(d1_i, d2_i), ...)."""
    d1_dict = dict(d1)
    result = {}
    for gen, exp in d2:
        if gen in d1_dict:
            result[gen] = max(exp, d1_dict[gen])
    return tuple(sorted(result.items()))


def max_dtuple(d1, d2):
    """return (max(d1_i, d2_i), ...)."""
    result = dict(d1)
    for gen, exp in d2:
        result[gen] = min(exp, result[gen]) if gen in result else exp
    return tuple(sorted(result.items()))


def gcd_nonzero_dtuple(d1, d2):
    """return (min(d1_i, d2_i), ...)."""
    d1_dict = dict(d1)
    for gen, exp in d2:
        if gen in d1_dict:
            return True
    return False
