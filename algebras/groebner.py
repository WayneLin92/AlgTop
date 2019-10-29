"""Algebras based on Groebner basis."""
# todo: create class AugModuleMod2
# Todo: construct subalgebra
import copy, itertools, operator, heapq, pickle
from typing import Tuple, List, Dict, Set,  Type, Iterable
from . import BaseAlgebras as BA, linalg, mymath


class GbAlgMod2(BA.AlgebraMod2):
    """A factory for algebras using Groebner basis.

    `GbAlgMod2.new_alg()` creates a new algebra with
    its own generators and relations.
    """

    _gen_names = None  # type: List[str]
    _gen_degs = None  # type: list
    _unit_deg = None
    _rels = None  # type: Dict[tuple, set]
    _rel_gen_leads = None  # type: Set[tuple]
    _key = None  # key function: mon -> value
    auto_simplify = None  # type: bool
    _name_index = 0

    @staticmethod
    def new_alg(*, unit_deg=None, key=None) -> "Type[GbAlgMod2]":
        """Return a dynamically created subclass of GbAlgMod2."""
        cls = GbAlgMod2
        class_name = f"GbAlgMod2_{cls._name_index}"
        cls._name_index += 1
        dct = {'_gen_names': [], '_gen_degs': [], '_unit_deg': unit_deg or 0,
               '_rels': {}, '_rel_gen_leads': set(), '_key': key, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbAlgMod2]":
        """Return a copy of current algebra."""
        class_name = f"GbAlgMod2_{GbAlgMod2._name_index}"
        GbAlgMod2._name_index += 1
        dct = {'_gen_names': cls._gen_names.copy(), '_gen_degs': cls._gen_degs.copy(),
               '_unit_deg': cls._unit_deg, '_rels': copy.deepcopy(cls._rels),
               '_rel_gen_leads': cls._rel_gen_leads.copy(), 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbAlgMod2,), dct)

    @classmethod
    def save(cls, filename):
        """Save to a pickle file."""
        with open(filename, 'wb') as file:
            pickle.dump([cls._gen_names, cls._gen_degs, cls._unit_deg,
                         cls._rels, cls._rel_gen_leads, cls._key], file)

    @staticmethod
    def load_alg(filename):
        """Create an algebra from a pickle file."""
        with open(filename, 'rb') as file:
            init_list = pickle.load(file)
            cls = GbAlgMod2
            class_name = f"GbAlgMod2_{cls._name_index}"
            cls._name_index += 1
            dct = {'_gen_names': init_list[0], '_gen_degs': init_list[1], '_unit_deg': init_list[2],
                   '_rels': init_list[3], '_rel_gen_leads': init_list[4], '_key': init_list[5], 'auto_simplify': True}
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
            return "".join(mymath.tex_pow(b, e) for b, e in zip(cls._gen_names, mon) if e)
        else:
            return "1"

    @classmethod
    def deg_mon(cls, mon: tuple):
        return sum(map(operator.mul, mon, cls._gen_degs), cls._unit_deg)

    def deg(self):
        """Require `self` to be homogeneous."""
        for m in self.data:
            return self.deg_mon(m)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, name: str, deg):
        """Add a new generator and return it."""
        cls._gen_names.append(name)
        cls._gen_degs.append(deg)
        m = (0,) * (len(cls._gen_names) - 1) + (1,)
        return cls(m).simplify()

    @classmethod
    def add_gens(cls, names_degs):
        """Add generators. names_degs is a list of tuples (name, deg)."""
        for nd in names_degs:
            cls._gen_names.append(nd[0])
            cls._gen_degs.append(nd[1])

    @classmethod
    def remove_gen(cls, name: str):
        """If `name=0` in the algebra with relations simplified, call this function to remove this generator."""
        i = cls._gen_names.index(name)
        m_k = (0,) * i + (1,)
        del cls._rels[m_k]

        def f(_m):
            return _m[:i] + _m[i+1:]
        cls._rels = {f(_m): {f(m1) for m1 in _v} for _m, _v in cls._rels.items()}
        cls._rel_gen_leads = set(map(f, cls._rel_gen_leads))
        cls._gen_names = f(cls._gen_names)
        cls._gen_degs = f(cls._gen_degs)

    @classmethod
    def rename_gen(cls, old_name, new_name):
        """Rename a generator."""
        i = cls._gen_names.index(old_name)
        cls._gen_names[i] = new_name

    @classmethod
    def reorder_gens(cls, index_map, key=None):
        """Reorganize the relations by a new ordering of generators and a new key function.
        The new i'th generator is the old `index_map[i]`'th generator."""
        num_gens = len(cls._gen_names)
        assert num_gens == len(index_map)

        def f(m):
            n = len(m)
            m1 = tuple(m[index_map[i]] if index_map[i] < n else 0 for i in range(num_gens))
            return mymath.rstrip_tuple(m1)
        rel_generators = cls.get_rel_gens()
        cls._key = key
        rel_generators = [{f(m) for m in rel} for rel in rel_generators]
        cls._gen_names = [cls._gen_names[index_map[i]] for i in range(num_gens)]
        cls._gen_degs = [cls._gen_degs[index_map[i]] for i in range(num_gens)]
        cls._rels = {}
        cls._rel_gen_leads = set()
        for rel in rel_generators:
            cls.add_rel_data(rel)

    @classmethod
    def add_rel(cls, rel: "GbAlgMod2"):
        """Add a relation."""
        if not rel.is_homo():
            raise ValueError(f'relation {rel} not homogeneous!')
        cls.add_rel_data(rel.data)

    @classmethod
    def add_rels(cls, rels: Iterable["GbAlgMod2"]):
        """Add a relation."""
        for rel in rels:
            if not rel.is_homo():
                raise ValueError(f'relation {rel} not homogeneous!')
        cls.add_rels_data(rel.data for rel in rels)

    @classmethod
    def simplify_rels(cls):
        """Simplify `cls._rels`."""
        for m in cls._rels:
            cls._rels[m] = cls.simplify_data(cls._rels[m])

    def simplify(self):
        """Simplify self by relations."""
        self.data = self.simplify_data(self.data)
        return self

    @classmethod
    def reduce_frob(cls):
        """Reduce the algebra by the square root of the ideal of relations."""
        # todo: fix this
        rels, cls._rels = cls._rels, {}
        n_gen = len(cls._gen_names)

        #
        import string
        cls._gen_names += [string.capwords(m) for m in cls._gen_names]
        cls._gen_degs += [2 * d for d in cls._gen_degs]

        for i in range(n_gen):
            y = (0,) * (n_gen + i) + (1,)
            x2 = (0,) * i + (2,)
            cls.add_rel_data({x2, y})
        # print("number of relations:", len(rels))
        for m, v in sorted(rels.items(), reverse=True):
            m_del = []
            for m1 in sorted(cls._rels, reverse=True):
                if m1 > m:
                    m_del.append(m1)
                else:
                    for m2 in m_del:
                        del cls._rels[m2]
                    break
            # print("rel:", cls(m), '=', cls(v), len(cls._rels))
            cls.add_rel(v | {m})
        m_del = []
        zero = (0,) * n_gen
        for m in sorted(cls._rels):
            if len(m) <= n_gen or m[:n_gen] != zero:
                m_del.append(m)

        for m in m_del:
            del cls._rels[m]
        cls._rels = dict((m[n_gen:], {_m[n_gen:] for _m in v})
                         for m, v in cls._rels.items())
        return rels != cls._rels

    @classmethod
    def reduce(cls):
        """Reduce all nilpotent elements."""
        while cls.reduce_frob():
            pass

    # getters --------------------------
    @classmethod
    def deg_data(cls, data: set):
        """Return the degree of `data`."""
        for m in data:
            return cls.deg_mon(m)

    @classmethod
    def get_lead(cls, data):
        """Return the leading term of `data`."""
        return max(data, key=cls._key) if cls._key else max(data)

    @classmethod
    def get_num_gens(cls):
        """Return the number of generators."""
        return len(cls._gen_names)

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = cls._gen_names.index(k)
        m = (0,) * i + (1,)
        return cls(m).simplify() if cls.auto_simplify else cls(m)

    @classmethod
    def get_rel_gens(cls, ideal=None):
        """Return a minimal generating set of `cls._rels` or `ideal`, ordered by degree."""
        rel_gens = []
        if ideal is None:
            for m in sorted(cls._rel_gen_leads, key=cls.deg_mon):
                rel = cls._rels[m] | {m}
                rel_gens.append(rel)
            return rel_gens
        else:
            rels_backup = cls._rels.copy()
            try:
                def deg_data(data):
                    for _m in data:
                        return cls(_m).deg()
                for rel in sorted(ideal, key=deg_data):
                    rel = cls.simplify_data(rel)
                    if rel:
                        rel_gens.append(rel)
                        cls.add_rel_data(rel)
            finally:
                cls._rels = rels_backup
            return rel_gens

    @classmethod
    def basis_mons_max(cls, deg_max):
        """Return an list of basis."""
        result = [((), 0)]
        for k in range(len(cls._gen_degs)):
            length = len(result)
            for i in range(length):
                m, d = result[i]
                for e in range(1, (deg_max - d) // cls._gen_degs[k] + 1):
                    m1 = m + (0,) * (k - len(m)) + (e,)
                    if any(map(mymath.le_tuple, cls._rels, itertools.repeat(m1))):
                        break
                    else:
                        result.append((m1, d + e * cls._gen_degs[k]))
        return result

    @classmethod
    def basis_max(cls, deg_max):
        return (cls(m) for m, d in cls.basis_mons_max(deg_max))

    @classmethod
    def is_reducible(cls, mon):
        """Determine if mon is reducible by `cls._rels`."""
        return any(mymath.le_tuple(m, mon) for m in cls._rels)

    def evaluation(self, image_gens):
        """Return f(self) where f is an algebraic map determined by `image_gens`."""
        assert len(image_gens) > 0
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
    def print_tex(cls, show_gb=True):
        """For pdflatex."""
        print(f"Generators: ${', '.join(cls._gen_names)}$.\\\\")
        print(f"Degrees: ${', '.join(map(str, cls._gen_degs))}$")
        print("Relations:\\\\")
        if show_gb:
            for m in cls._rels:
                if m in cls._rel_gen_leads:
                    print(f"$\\triangle\\hspace{{4pt}}{cls(m)} = {cls(cls._rels[m])}$\\\\")
                else:
                    print(f"${cls(m)} = {cls(cls._rels[m])}$\\\\")
        else:
            for m in cls._rel_gen_leads:
                print(f"${cls(m)} = {cls(cls._rels[m])}$\\\\")

    @classmethod
    def display_alg(cls, show_gb=True):
        """For Jupyter notebook."""
        from IPython.display import Markdown
        tr1 = "<th>Generators</th>"
        for name in cls._gen_names:
            tr1 += f"<td>${name}$</td>"
        tr1 = f"<tr>{tr1}</tr>"

        tr2 = "<th>Degrees</th>"
        for deg in cls._gen_degs:
            tr2 += f"<td>${deg}$</td>"
        tr2 = f"<tr>{tr2}</tr>"

        if show_gb:
            tr3 = "<th>Groebner basis</th>"
            td = "\\begin{align*}"
            for m in sorted(cls._rels, key=cls.deg_mon):
                if m in cls._rel_gen_leads:
                    td += "\\bullet\\hspace{4pt}"
                else:
                    td += "\\square\\hspace{4pt}"
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\\n"
            td += "\\end{align*}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{align*}"
            for m in sorted(cls._rel_gen_leads, key=cls.deg_mon):
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\\n"
            td += "\\end{align*}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        result = "<table>" + tr1 + tr2 + "</table>" + "<table>" + tr3 + "</table>"
        return Markdown(result)

    @classmethod
    def display_ideal(cls, gb: Iterable[set]):
        from IPython.display import Markdown
        result = "\\begin{align*}"
        for data in gb:
            result += f"&{cls(data)}\\\\\n"
        result += "\\end{align*}"
        return Markdown(result)

    # algorithms --------------------------
    @classmethod
    def simplify_data(cls, data: set):
        """Return simplified `data`. `data` will remain unchanged."""
        state_backup = cls.auto_simplify
        cls.auto_simplify = False
        s = data.copy()
        result = set()
        leading_masks = tuple({i for i, e in enumerate(m) if e} for m in cls._rels)
        while s:
            mon = max(s, key=cls._key) if cls._key else max(s)
            s.remove(mon)
            mask_mon = {i for i, e in enumerate(mon) if e}
            for m, mask_m in zip(cls._rels, leading_masks):
                if mask_m <= mask_mon and mymath.le_tuple(m, mon):
                    q, r = mymath.div_mod_tuple(mon, m)
                    m_to_q = (cls(cls._rels[m]) ** q).data
                    s ^= {mymath.add_tuple(r, m1) for m1 in m_to_q}
                    break
            else:
                result.add(mon)
        cls.auto_simplify = state_backup
        return result

    @classmethod
    def add_rel_data(cls, rel: set):
        """Add a relation."""
        rel = cls.simplify_data(rel)
        if not rel:
            return
        m = max(rel, key=cls._key) if cls._key else max(rel)
        cls._rel_gen_leads.add(m)
        hq = [(cls.deg_mon(m), rel)]
        while hq:
            _, r = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                m = max(r, key=cls._key) if cls._key else max(r)
                redundant_leading_terms = []
                for m1, v1 in cls._rels.items():
                    if any(map(min, m, m1)):  # gcd > 0
                        if mymath.le_tuple(m, m1):
                            redundant_leading_terms.append(m1)
                        lcm = mymath.max_tuple(m, m1)
                        dif = mymath.sub_tuple(lcm, m)
                        dif1 = mymath.sub_tuple(lcm, m1)
                        new_rel = {mymath.add_tuple(_m, dif) for _m in r}
                        v1dif1 = {mymath.add_tuple(_m, dif1) for _m in v1}
                        new_rel -= {lcm}
                        new_rel ^= v1dif1
                        heapq.heappush(hq, (cls.deg_mon(lcm), new_rel))
                for m_redundant in redundant_leading_terms:
                    del cls._rels[m_redundant]
                    if m_redundant in cls._rel_gen_leads:
                        cls._rel_gen_leads.remove(m_redundant)
                cls._rels[m] = r - {m}

    @classmethod
    def add_rels_data(cls, rels: Iterable[set]):
        """Add relations."""
        hq = [(cls.deg_data(rel), rel, True) for rel in rels]
        heapq.heapify(hq)
        while hq:
            # print(len(hq))  # ###############
            _, r, is_rel_gen = heapq.heappop(hq)
            r = cls.simplify_data(r)
            if r:
                if is_rel_gen:
                    m = max(r, key=cls._key) if cls._key else max(r)
                    cls._rel_gen_leads.add(m)
                m = max(r, key=cls._key) if cls._key else max(r)
                redundant_leading_terms = []
                for m1, v1 in cls._rels.items():
                    if any(map(min, m, m1)):  # gcd > 0
                        if mymath.le_tuple(m, m1):
                            redundant_leading_terms.append(m1)
                        lcm = mymath.max_tuple(m, m1)
                        dif = mymath.sub_tuple(lcm, m)
                        dif1 = mymath.sub_tuple(lcm, m1)
                        new_rel = {mymath.add_tuple(_m, dif) for _m in r}
                        v1dif1 = {mymath.add_tuple(_m, dif1) for _m in v1}
                        new_rel -= {lcm}
                        new_rel ^= v1dif1
                        heapq.heappush(hq, (cls.deg_mon(lcm), new_rel, False))
                for m_redundant in redundant_leading_terms:
                    del cls._rels[m_redundant]
                    if m_redundant in cls._rel_gen_leads:
                        cls._rel_gen_leads.remove(m_redundant)
                cls._rels[m] = r - {m}

    @classmethod
    def ann(cls, x):
        """Return the groebner basis for the ideal {y | xy=0}."""
        rels_backup = cls._rels.copy()
        rel_gen_leads_backup = cls._rel_gen_leads.copy()
        key_backup = cls._key
        d = x.deg()
        try:
            X = cls.add_gen('X_{ann}', d)
            num_gen = len(cls._gen_names)
            if key_backup is None:
                cls._key = lambda _m: (-_m[-1] if len(_m) == num_gen else 0, _m)
            else:
                cls._key = lambda _m: (-_m[-1] if len(_m) == num_gen else 0, key_backup(_m))
            cls.add_rel(X + x)
            rels = cls._rels
        finally:
            cls._gen_names.pop()
            cls._gen_degs.pop()
            cls._rels = rels_backup
            cls._rel_gen_leads = rel_gen_leads_backup
            cls._key = key_backup
        result = []
        for m in rels:
            if len(m) == num_gen:
                result_i = cls.zero()
                for m1 in itertools.chain((m,), rels[m]):
                    result_i += cls(mymath.rstrip_tuple(m1[:-1])) * (x ** (m1[-1] - 1))
                result.append(result_i.data)
        return result

    @classmethod
    def ann_seq(cls, ele_names: List[Tuple["GbAlgMod2", str]]):
        """Display relations among elements: $\\sum a_ix_i=0$."""
        from IPython.display import Markdown
        rels_backup = cls._rels.copy()
        rel_gen_leads_backup = cls._rel_gen_leads.copy()
        key_backup = cls._key
        try:
            num_gen = len(cls._gen_names)
            if key_backup is None:
                cls._key = lambda _m: ([-i for i in _m[num_gen:]] + [0] * (len(ele_names) - len(_m[num_gen:])), _m)
            else:
                cls._key = lambda _m: ([-i for i in _m[num_gen:]] + [0] * (len(ele_names) - len(_m[num_gen:])),
                                       key_backup(_m))
            for ele, name in ele_names:
                x = cls.add_gen(name, ele.deg())
                cls.add_rel(x + ele)
            result = ""
            for m in cls._rels:
                if len(m) > num_gen:
                    result += f"* ${cls({m} | cls._rels[m])} = 0$\n"
        finally:
            cls._gen_names.pop()
            cls._gen_degs.pop()
            cls._rels = rels_backup
            cls._rel_gen_leads = rel_gen_leads_backup
            cls._key = key_backup
        return Markdown(result)

    @classmethod
    def subalgebra(cls, ele_names: List[Tuple["GbAlgMod2", str]], *, key=None):
        """Return the subalgebra generated by `ele_names`."""
        assert cls._key is None
        num_gens = len(cls._gen_names)
        A = cls.copy_alg()
        if key:
            def key1(_m):
                return _m[:num_gens], key(_m[num_gens:])
        else:
            def key1(_m):
                return _m
        A._key = key1
        for ele, name in ele_names:
            x = A.add_gen(name, ele.deg())
            A.add_rel(x + ele)
        A._gen_names = A._gen_names[num_gens:]
        A._gen_degs = A._gen_degs[num_gens:]
        A._key = key
        A._rel_gen_leads = set()
        rels, A._rels = A._rels, {}
        for m in rels:
            n = 0
            while n < len(m) and m[n] == 0:
                n += 1
            if n >= num_gens:
                A.add_rel_data({m[num_gens:]} | {_m[num_gens:] for _m in rels[m]})
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
        dct = {'_gen_names': [], '_gen_degs': [], '_gen_diff': [], '_unit_deg': unit_deg or 0,
               '_rels': {}, '_rel_gen_leads': set(), '_key': key, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbDga]":
        """Return a copy of current algebra."""
        class_name = f"GbDGA_{GbDga._name_index}"
        GbDga._name_index += 1
        dct = {'_gen_names': cls._gen_names.copy(), '_gen_degs': cls._gen_degs.copy(),
               '_gen_diff': cls._gen_diff, '_unit_deg': cls._unit_deg,
               '_rels': copy.deepcopy(cls._rels), 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbDga,), dct)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, k: str, deg, diff=None):
        """Add a new generator and return it."""
        cls._gen_names.append(k)
        cls._gen_degs.append(deg)
        if diff is None:
            diff = set()
        elif type(diff) is not set:
            diff = diff.data
        cls._gen_diff.append(diff)
        m = (0,) * (len(cls._gen_names) - 1) + (1,)
        return cls(m).simplify()

    def diff(self):
        """Return the coboundary of the cochain."""
        result = set()
        for m in self.data:
            for i in range(len(m)):
                if m[i] % 2:
                    m1 = mymath.rstrip_tuple(m[:i] + (m[i] - 1,) + m[i+1:])
                    m1_by_dg_i = {mymath.add_tuple(m1, _m) for _m in self._gen_diff[i]}
                    result ^= m1_by_dg_i
        return type(self)(result).simplify()

    # getters ----------------------------
    @classmethod
    def homology(cls, deg_max) -> Tuple[Type[GbAlgMod2], list]:
        """Compute HA. Return (HA, list of representing cycles)."""
        map_diff = linalg.GradedLinearMapKMod2()
        for r in cls.basis_max(deg_max):
            map_diff.add_map(r, r.diff())
        Z = [map_diff.kernel(d) for d in range(deg_max + 1)]
        B = [map_diff.image(d) for d in range(deg_max + 1)]
        H = [Z[d] / B[d] for d in range(deg_max + 1)]

        R = GbAlgMod2.new_alg()
        R_basis_mons = [((), 0)]
        map_alg = linalg.GradedLinearMapKMod2()
        image_gens = []
        index = 1
        for d in range(1, deg_max + 1):
            for x in (H[d] / map_alg.image(d)).basis(cls):
                R.add_gen(mymath.tex_sub('x', index), d)
                index += 1
                image_gens.append(x)

                length = len(R_basis_mons)
                for i in range(length):
                    m1, d1 = R_basis_mons[i]
                    for e in range(1, (deg_max - d1) // d + 1):
                        m2 = m1 + (0,) * (len(R._gen_degs) - len(m1) - 1) + (e,)
                        d2 = d1 + e * d
                        R_basis_mons.append((m2, d2))
                        r2 = R(m2)
                        fr2 = B[d2].res(r2.evaluation(image_gens))
                        map_alg.add_map(r2, fr2)
                        if map_alg.kernel(d2):
                            for rel in map_alg.kernel(d2).basis(R):
                                R.add_rel(rel)
                            map_alg.kernel(d2).clear()
        return R, image_gens

    @classmethod
    def resolution(cls, deg_max) -> Type["GbDga"]:
        """Compute Tor_A(k, k)."""
        R = cls.copy_alg()
        R_basis_mons = R.basis_mons_max(deg_max)
        map_diff = linalg.GradedLinearMapKMod2()
        for m, d in R_basis_mons:
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
                    m1, d1 = R_basis_mons[i]
                    if d1 + d <= deg_max:
                        m2 = m1 + (0,) * (len(R._gen_degs) - len(m1) - 1) + (1,)
                        d2 = d1 + d
                        R_basis_mons.append((m2, d2))
                        r2 = R(m2)
                        map_diff.add_map(r2, r2.diff())
        return R

    @classmethod
    def print_tex(cls, show_gb=False):
        """Print the cls in latex."""
        super().print_tex(show_gb)
        print("Differentials:\\\\")
        for g, dg in zip(cls._gen_names, cls._gen_diff):
            print(f"$d({g})={cls(dg)}$\\\\")

    @classmethod
    def display_alg(cls, show_gb=False):
        from IPython.display import Markdown
        td_style = 'style="text-align:left;"'
        tr1 = '<th>Generators</th>'
        for name in cls._gen_names:
            tr1 += f'<td {td_style}>${name}$</td>'
        tr1 = f'<tr>{tr1}</tr>\n'

        tr2 = '<th>Generators</th>'
        for deg in cls._gen_degs:
            tr2 += f'<td {td_style}>${deg}$</td>'
        tr2 = f'<tr>{tr2}</tr>\n'

        if show_gb:
            tr3 = "<th>Groebner basis</th>"
            td = "\\begin{align*}"
            for m in sorted(cls._rels, key=cls.deg_mon):
                if m in cls._rel_gen_leads:
                    td += "\\bullet\\hspace{4pt}"
                else:
                    td += "\\square\\hspace{4pt}"
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\\n"
            td += "\\end{align*}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{align*}"
            for m in sorted(cls._rel_gen_leads, key=cls.deg_mon):
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\\n"
            td += "\\end{align*}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        tr4 = '<th>Differentials</th>'
        td = '\\begin{align*}'
        for g, dg in zip(cls._gen_names, cls._gen_diff):
            td += f'd({g}) &= {cls(dg)}\\\\\n'
        td += '\\end{align*}'
        td = f'<td>{td}</td>'
        tr4 += td
        tr4 = f'<tr>{tr4}</tr>\n'

        result = '<table>\n' + tr1 + tr2 + '</table>' +\
                 '<table>' + tr3 + tr4 + '</table>'
        return Markdown(result)
