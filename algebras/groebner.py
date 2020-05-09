"""Algebras based on Groebner basis."""
# TODO: get rid of unit_deg
import copy
import heapq
from itertools import chain, repeat, combinations, groupby
import operator
import pickle
from typing import Tuple, List, Dict, Set, Type, Iterable, Callable, Any
from algebras import BaseAlgebras as BA, linalg, mymath


class GbAlgMod2(BA.AlgebraMod2):
    """A factory for algebras using Groebner basis.

    `GbAlgMod2.new_alg()` creates a new algebra with
    its own generators and relations.
    """

    gen_names = None  # type: List[str]
    gen_degs = None  # type: list
    _unit_deg = None
    rels = None  # type: Dict[tuple, Set[Tuple[int]]]
    _rels_gen_leads = None  # type: Set[tuple]
    _rels_cache = None  # type: List[Tuple[int, bool, set]]
    key = None  # type: Callable[[tuple], Any]
    auto_simplify = None  # type: bool
    _attributes = ["gen_names", "gen_degs", "_unit_deg", "rels",
                   "_rels_gen_leads", "_rels_cache", "key", "auto_simplify"]
    _name_index = 0

    @staticmethod
    def new_alg(*, unit_deg=None, key=None) -> "Type[GbAlgMod2]":
        """Return a dynamically created subclass of GbAlgMod2.

        When key=None, use reversed lexicographical ordering by default."""
        cls = GbAlgMod2
        class_name = f"GbAlgMod2_{cls._name_index}"
        cls._name_index += 1
        dct = {'gen_names': [], 'gen_degs': [], '_unit_deg': unit_deg or 0,
               'rels': {}, '_rels_gen_leads': set(), '_rels_cache': [],
               'key': key, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbAlgMod2]":
        """Return a copy of current algebra."""
        class_name = f"GbAlgMod2_{GbAlgMod2._name_index}"
        GbAlgMod2._name_index += 1
        dct = {'gen_names': cls.gen_names.copy(), 'gen_degs': cls.gen_degs.copy(),
               '_unit_deg': cls._unit_deg, 'rels': copy.deepcopy(cls.rels),
               '_rels_gen_leads': cls._rels_gen_leads.copy(),
               '_rels_cache': copy.deepcopy(cls._rels_cache),
               'key': cls.key, 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbAlgMod2,), dct)

    @classmethod
    def save_alg(cls, filename):
        """Save to a pickle file."""
        with open(filename, 'wb') as file:
            pickle.dump([getattr(cls, attr) for attr in cls._attributes], file)

    @staticmethod
    def load_alg(filename):  # TODO: change this
        """Create an algebra from a pickle file."""
        with open(filename, 'rb') as file:
            init_list = pickle.load(file)
            cls = GbAlgMod2
            class_name = f"GbAlgMod2_{cls._name_index}"
            cls._name_index += 1
            attributes = ["gen_names", "gen_degs", "_unit_deg", "rels",
                          "_rels_gen_leads", "key"]
            dct = {attr: init_v for attr, init_v in zip(cls._attributes, init_list)}
            dct["_rels_cache"] = []
            dct["auto_simplify"] = True
            # noinspection PyTypeChecker
            return type(class_name, (cls,), dct)

    # ----- AlgebraMod2 -------------
    @classmethod
    def mul_mons(cls, mon1: tuple, mon2: tuple):
        m = mymath.add_tuple(mon1, mon2)
        return cls.simplify_data({m}) if cls.auto_simplify else m

    @classmethod
    def str_mon(cls, mon: tuple):
        if mon:
            return "".join(mymath.tex_pow(b, e) for b, e in zip(cls.gen_names, mon) if e)
        else:
            return "1"

    @classmethod
    def repr_mon(cls, mon: tuple, clsname):
        if mon:
            return " * ".join(f"{clsname}.gen(\"{b}\") ** {e}" if e > 1 else f"{clsname}.gen(\"{b}\")"
                              for b, e in zip(cls.gen_names, mon) if e)
        else:
            return f"{clsname}.unit()"

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(operator.mul, mon, cls.gen_degs), cls._unit_deg)

    def deg(self):
        """Require `self` to be homogeneous."""
        for m in self.data:
            return self.deg_mon(m)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, name: str, deg):
        """Add a new generator and return it."""
        cls.gen_names.append(name)
        cls.gen_degs.append(deg)
        m = (0,) * (len(cls.gen_names) - 1) + (1,)
        return cls(m).simplify() if cls.auto_simplify else cls(m)

    @classmethod
    def add_gens(cls, names_degs):
        """Add generators. names_degs is a list of tuples (name, deg)."""
        for nd in names_degs:
            cls.gen_names.append(nd[0])
            cls.gen_degs.append(nd[1])

    @classmethod
    def remove_gen(cls, name: str):
        """Remove a generator.

        If the generator `name` equals zero in the algebra whose relations are simplified,
        call this function to remove this generator. Use only when _rels_cache is empty."""
        i = cls.gen_names.index(name)
        m_k = (0,) * i + (1,)
        assert not cls.rels[m_k]
        del cls.rels[m_k]

        def f(_m):
            return _m[:i] + _m[i+1:]
        cls.rels = {f(_m): {f(m1) for m1 in _v} for _m, _v in cls.rels.items()}
        cls._rels_gen_leads = {fm for _m in cls._rels_gen_leads if sum(fm := f(_m))}
        cls.gen_names = f(cls.gen_names)
        cls.gen_degs = f(cls.gen_degs)

    @classmethod
    def rename_gen(cls, old_name, new_name):
        """Rename a generator."""
        i = cls.gen_names.index(old_name)
        cls.gen_names[i] = new_name

    @classmethod
    def reorder_gens(cls, index_map=None, key=None):  # TODO: create a new alg instead
        """Reorganize the relations by a new ordering of generators and a new key function.
        The new i'th generator is the old `index_map[i]`'th generator."""
        num_gens = len(cls.gen_names)
        rel_generators = cls.get_rel_gens()
        cls.key = key
        if index_map:
            def f(m):
                n = len(m)
                m1 = tuple(m[index_map[i]] if index_map[i] < n else 0 for i in range(num_gens))
                return mymath.rstrip_tuple(m1)
            assert num_gens == len(index_map)
            rel_generators = [{f(m) for m in rel} for rel in rel_generators]
            cls.gen_names = [cls.gen_names[index_map[i]] for i in range(num_gens)]
            cls.gen_degs = [cls.gen_degs[index_map[i]] for i in range(num_gens)]
        cls.rels = {}
        cls._rels_gen_leads = set()
        cls.add_rels_data(rel_generators, clear_cache=True)

    @classmethod
    def add_rel_data(cls, rel: set, clear_cache=False):
        """Add relations."""
        if rel:
            deg = cls.deg_data(rel)
            heapq.heappush(cls._rels_cache, (deg, True, rel))
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
        for i, item in enumerate(cls._rels_cache):
            cls._rels_cache[i] = (item[0], item[1], cls.simplify_data(item[2]))

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    # getters --------------------------
    @classmethod
    def deg_data(cls, data: set):
        """Return the degree of `data`."""
        for m in data:
            return cls.deg_mon(m)

    @classmethod
    def get_lead(cls, data):
        """Return the leading term of `data`."""
        return max(data, key=cls.key) if cls.key else min(data)

    @classmethod
    def get_num_gens(cls):
        """Return the number of generators."""
        return len(cls.gen_names)

    @classmethod
    def get_cache_size(cls):
        return len(cls._rels_cache)

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = cls.gen_names.index(k)
        m = (0,) * i + (1,)
        return cls(m).simplify() if cls.auto_simplify else cls(m)

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
    def basis_mons_max(cls, deg_max):
        """Return a list of basis (mon, deg)."""
        result = [(cls._unit_deg, ())]
        leadings = sorted(cls.rels, key=len)
        leadings = {k - 1: list(g) for k, g in groupby(leadings, key=len)}
        for k in range(len(cls.gen_degs)):
            for i in range(len(result)):
                d, m = result[i]
                for e in range(1, (deg_max - d) // cls.gen_degs[k] + 1):
                    m1 = m + (0,) * (k - len(m)) + (e,)
                    if k in leadings and any(map(mymath.le_tuple, leadings[k], repeat(m1))):
                        break
                    else:
                        result.append((d + e * cls.gen_degs[k], m1))
        return result

    @classmethod
    def basis_max(cls, deg_max):
        return ((d, cls(m)) for d, m in cls.basis_mons_max(deg_max))

    @classmethod
    def is_reducible(cls, mon):
        """Determine if mon is reducible by `cls.rels`."""
        return any(mymath.le_tuple(m, mon) for m in cls.rels)

    def evaluation(self, image_gens):
        """Return f(self) where f is an algebraic map determined by `image_gens`."""
        assert len(image_gens) == len(self.gen_names)
        R = type(image_gens[0])
        zero = R.zero() if issubclass(R, BA.Algebra) else 0
        unit = R.unit() if issubclass(R, BA.Algebra) else 1
        result = zero
        for m in self.data:
            fm = unit
            for fg, e in zip(image_gens, m):
                fm *= (fg ** e)
            result += fm
        return result

    @classmethod
    def print_latex_alg(cls, show_gb=True):
        """For latex."""
        print(f"Generators: ${', '.join(cls.gen_names)}$.\\\\")
        print(f"Degrees: ${', '.join(map(str, cls.gen_degs))}$")
        print("Relations:\\\\")
        if show_gb:
            for m in cls.rels:
                if m in cls._rels_gen_leads:
                    print(f"${cls(m)} = {cls(cls.rels[m])}$\n")
                else:
                    print(f"$\\bullet {cls(m)} = {cls(cls.rels[m])}$\n")
        else:
            for m in cls._rels_gen_leads:
                print(f"${cls(m)} = {cls(cls.rels[m])}$\\\\")

    @classmethod
    def markdown_alg(cls, show_gb=True):
        """For Jupyter notebook."""
        from IPython.display import Markdown
        tr1 = "<th>Generators</th>"
        for name in cls.gen_names:
            tr1 += f"<td>${name}$</td>"
        tr1 = f"<tr>{tr1}</tr>"

        tr2 = "<th>Degrees</th>"
        for deg in cls.gen_degs:
            tr2 += f"<td>${deg}$</td>"
        tr2 = f"<tr>{tr2}</tr>"

        if show_gb:
            tr3 = "<th>Groebner basis</th>"
            td = "\\begin{align*}"
            for m in sorted(cls.rels, key=cls.deg_mon):
                if m in cls._rels_gen_leads:
                    td += "\\bullet\\hspace{4pt}"
                else:
                    td += "\\square\\hspace{4pt}"
                td += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            td += "\\end{align*}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{align*}\n"
            for m in sorted(cls._rels_gen_leads, key=cls.deg_mon):
                td += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            td += "\\end{align*}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        result = "<table>" + tr1 + tr2 + "</table>" + "<table>" + tr3 + "</table>"
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
            s = "+".join(f"{name}{mymath.tex_parenthesis(coeff)}" for coeff, name, d in a)
            result += f"& {s}=0\\\\\n"
        result += "\\end{align*}"
        return Latex(result)

    @staticmethod
    def repr_annilators(clsname, annilators:  List[List[Tuple["GbAlgMod2", str, int]]]):
        result = ""
        for a in annilators:
            s = " + ".join(f"{clsname}.gen(\"{name}\") * {mymath.tex_parenthesis(coeff.repr_(clsname))}"
                           for coeff, name, d in a)
            result += f"{s}\n"
        return result

    @classmethod
    def to_dga(cls) -> "Type[GbDga]":
        class_name = f"GbDGA_{GbDga._name_index}"
        GbDga._name_index += 1
        dct = {'gen_names': cls.gen_names.copy(), 'gen_degs': cls.gen_degs.copy(),
               '_gen_diff': [None] * len(cls.gen_names), '_unit_deg': cls._unit_deg,
               'rels': copy.deepcopy(cls.rels), 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbDga,), dct)

    # algorithms --------------------------
    @classmethod
    def simplify_data(cls, data: set):
        """Return simplified `data`. `data` will remain unchanged."""
        state_backup = cls.auto_simplify
        cls.auto_simplify = False
        s = data.copy()
        result = set()
        leading_masks = tuple({i for i, e in enumerate(m) if e} for m in cls.rels)
        while s:
            mon = cls.get_lead(s)
            s.remove(mon)
            mask_mon = {i for i, e in enumerate(mon) if e}
            for m, mask_m in zip(cls.rels, leading_masks):
                if mask_m <= mask_mon and mymath.le_tuple(m, mon):
                    q, r = mymath.div_mod_tuple(mon, m)
                    m_to_q = (cls(cls.rels[m]) ** q).data
                    s ^= {mymath.add_tuple(r, m1) for m1 in m_to_q}
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
            deg, is_rel_gen, r = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                m = cls.get_lead(r)
                if is_rel_gen:
                    cls._rels_gen_leads.add(m)
                redundant_leading_terms = []
                for m1, v1 in cls.rels.items():
                    if any(map(min, m, m1)):  # gcd > 0
                        lcm = mymath.max_tuple(m, m1)
                        dif = mymath.sub_tuple(lcm, m)
                        dif1 = mymath.sub_tuple(lcm, m1)
                        new_rel = {mymath.add_tuple(_m, dif) for _m in r}
                        v1dif1 = {mymath.add_tuple(_m, dif1) for _m in v1}
                        new_rel -= {lcm}
                        new_rel ^= v1dif1
                        if mymath.le_tuple(m, m1):
                            BA.Monitor.count("redundant_leading_terms")  #
                            redundant_leading_terms.append(m1)
                            is_lead = m1 in cls._rels_gen_leads
                            if new_rel:
                                heapq.heappush(hq, (cls.deg_mon(lcm), is_lead, new_rel))
                            if is_lead:
                                cls._rels_gen_leads.remove(m1)
                        elif new_rel:
                            heapq.heappush(hq, (cls.deg_mon(lcm), False, new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls.rels[m_redundant]
                cls.rels[m] = r - {m}

    @classmethod
    def get_vector_gens(cls, ideal: List[List[Tuple["GbAlgMod2", str, int]]], *, inplace=False):
        """Return the minimal generating set of `ideal`, which is a A-submodule of A^n.

        `ideal` is a list of [(ele, name, deg), ...].
        The names should not overlap with existing generator names of `cls`."""
        A = cls if inplace else cls.copy_alg()
        num_gen = len(A.gen_names)
        rels = []
        for index, v in enumerate(ideal):
            rel = A.zero()
            for ele, name, deg in v:
                name = f"v_{{{name}}}"
                x = A.gen(name) if name in A.gen_names else A.add_gen(name, deg)  #
                rel += x * ele
            if rel:
                rels.append((index, rel))
        num_gen1 = len(A.gen_names)
        rels_module = []
        for i in range(num_gen, num_gen1):
            for j in range(i, num_gen1):
                if i == j:
                    m = (0,) * i + (2,)
                else:
                    m = tuple(1 if k in (i, j) else 0 for k in range(len(A.gen_names)))
                rels_module.append(A(m))
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
        num_gen = len(cls.gen_names)
        num_ele = len(ele_names)
        if cls.key:
            A.key = lambda _m: ([-i for i in _m[num_gen:]] + [0] * num_ele, cls.key(_m))
        else:
            A.key = lambda _m: ([-i for i in _m[num_gen:]] + [0] * num_ele, [-i for i in _m])  #
        rels_new = []
        for ele, name in ele_names:
            x = A.add_gen(name, ele.deg())
            rels_new.append(x + ele)
        A.add_rels(rels_new, clear_cache=True)
        annilators = []
        for m in A.rels:
            if len(m) > num_gen:
                a = []
                for m1 in chain((m,), A.rels[m]):
                    name = A.gen_names[len(m1) - 1]
                    deg = A.gen_degs[len(m1) - 1]
                    m11 = mymath.rstrip_tuple(m1[:-1] + (m1[-1] - 1,))
                    a.append((m11, name, deg))
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
                return A.deg_mon(_m[num_gen:]), cls.key(_m)
        else:
            def key(_m):
                return A.deg_mon(_m[num_gen:]), [-i for i in _m]
        A.reorder_gens(key=key)
        annilators = [[(cls(A.simplify_data({_m})), name, deg) for _m, name, deg in a] for a in annilators]
        return cls.get_vector_gens(annilators)

    @classmethod
    def subalgebra(cls, ele_names: List[Tuple["GbAlgMod2", str]], *, key=None):
        """Return the subalgebra generated by `ele_names`."""
        num_gens = len(cls.gen_names)
        A = cls.copy_alg()

        def key1(_m):
            _m1, _m2 = _m[:num_gens], _m[num_gens:]
            return (cls.deg_mon(_m1), (cls.key(_m1) if cls.key else [-i for i in _m1]),
                    (key(_m2) if key else [-i for i in _m2]))
        A.key = key1
        for ele, name in ele_names:
            x = A.add_gen(name, ele.deg())
            A.add_rel(x + ele)
        A.add_rels_cache()
        A.gen_names = A.gen_names[num_gens:]
        A.gen_degs = A.gen_degs[num_gens:]
        A.key = key
        A._rels_gen_leads = set()
        rels, A.rels = A.rels, {}
        for m in rels:
            n = 0
            while n < len(m) and m[n] == 0:
                n += 1
            if n >= num_gens:
                rel_subalg = {m[num_gens:]} | {_m[num_gens:] for _m in rels[m]}
                A.add_rel_data(rel_subalg)
        A.add_rels_cache()
        return A


class GbDga(GbAlgMod2):
    """A factory for DGA over F_2."""
    _gen_diff = None  # type: List[set]

    @staticmethod
    def new_alg(*, unit_deg=None, key=None) -> "Type[GbDga]":
        """Return a dynamically created subclass of GbDga."""
        cls = GbDga
        class_name = f"GbDGA_{cls._name_index}"
        cls._name_index += 1
        dct = {'gen_names': [], 'gen_degs': [], '_gen_diff': [], '_unit_deg': unit_deg or 0,
               'rels': {}, '_rels_gen_leads': set(), 'key': key, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbDga]":
        """Return a copy of current algebra."""
        class_name = f"GbDGA_{GbDga._name_index}"
        GbDga._name_index += 1
        dct = {'gen_names': cls.gen_names.copy(), 'gen_degs': cls.gen_degs.copy(),
               '_gen_diff': cls._gen_diff, '_unit_deg': cls._unit_deg,
               'rels': copy.deepcopy(cls.rels), 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbDga,), dct)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, k: str, deg, diff=None):
        """Add a new generator and return it."""
        cls.gen_names.append(k)
        cls.gen_degs.append(deg)
        if diff is None:
            diff = set()
        elif type(diff) is not set:
            diff = diff.data
        cls._gen_diff.append(diff)
        m = (0,) * (len(cls.gen_names) - 1) + (1,)
        return cls(m).simplify()

    @classmethod
    def set_diff(cls, k: str, diff):
        index = cls.gen_names.index(k)
        cls._gen_diff[index] = diff.data

    def diff(self):
        """Return the boundary of the chain."""
        result = set()
        for m in self.data:
            for i in range(len(m)):
                if m[i] % 2:
                    m1 = mymath.rstrip_tuple(m[:i] + (m[i] - 1,) + m[i+1:])
                    m1_by_dg_i = {mymath.add_tuple(m1, _m) for _m in self._gen_diff[i]}
                    result ^= m1_by_dg_i
        return type(self)(result).simplify()

    def inv_diff(self):
        """Return the boundary of the chain."""
        pass

    # getters ----------------------------
    @classmethod
    def homology(cls, deg_max) -> Type[GbAlgMod2]:
        """Compute HA. Return (HA, list of representing cycles).

        HX._rel_cache is not cleared above deg_max."""
        map_diff = linalg.GradedLinearMapKMod2()
        for d, r in cls.basis_max(deg_max):
            # noinspection PyUnresolvedReferences
            map_diff.add_map(r, r.diff())
        Z = [map_diff.kernel(d) for d in range(deg_max + 1)]
        B = [map_diff.image(d) for d in range(deg_max + 1)]
        H = [Z[d] / B[d] for d in range(deg_max + 1)]

        R = GbAlgMod2.new_alg()
        R_basis_mons = [((), 0)]
        map_alg = linalg.GradedLinearMapKMod2()
        image_gens = []
        for d in range(1, deg_max + 1):
            for x in (H[d] / map_alg.image(d)).basis(cls):
                R.add_gen(f'[{x}]', d)
                image_gens.append(x)

                length = len(R_basis_mons)
                for i in range(length):
                    m1, d1 = R_basis_mons[i]
                    for e in range(1, (deg_max - d1) // d + 1):
                        m2 = m1 + (0,) * (len(R.gen_degs) - len(m1) - 1) + (e,)
                        d2 = d1 + e * d
                        R_basis_mons.append((m2, d2))
                        r2 = R(m2)
                        fr2 = B[d2].res(r2.evaluation(image_gens))
                        map_alg.add_map(r2, fr2)
                        if map_alg.kernel(d2):
                            R.add_rels_data(map_alg.kernel(d2).basis(set))
                            map_alg.kernel(d2).clear()

        return R

    @classmethod
    def resolution(cls, deg_max) -> Type["GbDga"]:
        """Compute Tor_A(k, k)."""
        R = cls.copy_alg()
        R_basis_mons = R.basis_mons_max(deg_max)
        map_diff = linalg.GradedLinearMapKMod2()
        for d, m in R_basis_mons:
            r = R(m)
            map_diff.add_map(r, r.diff())
        index = 1
        for d in range(1, deg_max + 1):
            h = map_diff.kernel(d) / map_diff.image(d)
            for x in h.basis(R):
                y = R.add_gen(mymath.tex_sub('y', index), d, x)
                R.add_rel(y * y)
                index += 1

                length = len(R_basis_mons)
                for i in range(length):
                    d1, m1 = R_basis_mons[i]
                    if d1 + d <= deg_max:
                        m2 = m1 + (0,) * (len(R.gen_degs) - len(m1) - 1) + (1,)
                        d2 = d1 + d
                        R_basis_mons.append((d2, m2))
                        r2 = R(m2)
                        map_diff.add_map(r2, r2.diff())
        return R

    @classmethod
    def print_latex_alg(cls, show_gb=False):
        """Print the cls in latex."""
        super().print_latex_alg(show_gb)
        print("Differentials:\\\\")
        for g, dg in zip(cls.gen_names, cls._gen_diff):
            print(f"$d({g})={cls(dg)}$\\\\")

    @classmethod
    def markdown_alg(cls, show_gb=False):
        from IPython.display import Markdown
        td_style = 'style="text-align:left;"'
        tr1 = '<th>Generators</th>'
        for name in cls.gen_names:
            tr1 += f'<td {td_style}>${name}$</td>'
        tr1 = f'<tr>{tr1}</tr>\n'

        tr2 = '<th>Generators</th>'
        for deg in cls.gen_degs:
            tr2 += f'<td {td_style}>${deg}$</td>'
        tr2 = f'<tr>{tr2}</tr>\n'

        if show_gb:
            tr3 = "<th>Groebner basis</th>"
            td = "\\begin{align*}"
            for m in sorted(cls.rels, key=cls.deg_mon):
                if m in cls._rels_gen_leads:
                    td += "\\bullet\\hspace{4pt}"
                else:
                    td += "\\square\\hspace{4pt}"
                td += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            td += "\\end{align*}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{align*}"
            for m in sorted(cls._rels_gen_leads, key=cls.deg_mon):
                td += f"{cls(m)} &= {cls(cls.rels[m])}\\\\\n"
            td += "\\end{align*}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        tr4 = '<th>Differentials</th>'
        td = '\\begin{align*}'
        for g, dg in zip(cls.gen_names, cls._gen_diff):
            td += f'd({g}) &= {cls(dg)}\\\\\n'
        td += '\\end{align*}'
        td = f'<td>{td}</td>'
        tr4 += td
        tr4 = f'<tr>{tr4}</tr>\n'

        result = '<table>\n' + tr1 + tr2 + '</table>' +\
                 '<table>' + tr3 + tr4 + '</table>'
        return Markdown(result)

# 691, 656, 737
