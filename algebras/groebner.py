"""Algebras based on Groebner basis."""
# todo: create class AugModuleMod2
# Todo: construct subalgebra
import copy, itertools, operator, heapq, pickle
from typing import Union, Set, Tuple, List, Dict, Type
from algebras import BaseAlgebras as BA, linalg, mymath


class GbAlgMod2(BA.AlgebraMod2):
    """A factory for algebras using Groebner basis.

    `GbAlgMod2.new_alg()` creates a new algebra with
    its own generators and relations.
    """

    _gen_names = None  # type: List[str]
    _gen_degs = None  # type: list
    _unit_deg = None
    _rels = None  # type: Dict[tuple, set]
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
               '_rels': {}, '_key': key, 'auto_simplify': True}
        # noinspection PyTypeChecker
        return type(class_name, (cls,), dct)

    @classmethod
    def copy_alg(cls) -> "Type[GbAlgMod2]":
        """Return a copy of current algebra."""
        class_name = f"GbAlgMod2_{GbAlgMod2._name_index}"
        GbAlgMod2._name_index += 1
        dct = {'_gen_names': cls._gen_names.copy(), '_gen_degs': cls._gen_degs.copy(),
               '_unit_deg': cls._unit_deg, '_rels': copy.deepcopy(cls._rels), 'auto_simplify': cls.auto_simplify}
        # noinspection PyTypeChecker
        return type(class_name, (GbAlgMod2,), dct)

    @classmethod
    def save(cls, filename):
        """Save to a pickle file."""
        with open(filename, 'wb') as file:
            pickle.dump([cls._gen_names, cls._gen_degs, cls._unit_deg, cls._rels, cls._key], file)

    @staticmethod
    def load_alg(filename):
        """Create an algebra from a pickle file."""
        with open(filename, 'rb') as file:
            init_list = pickle.load(file)
            cls = GbAlgMod2
            class_name = f"GbAlgMod2_{cls._name_index}"
            cls._name_index += 1
            dct = {'_gen_names': init_list[0], '_gen_degs': init_list[1], '_unit_deg': init_list[2],
                   '_rels': init_list[3], '_key': init_list[4], 'auto_simplify': True}
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
        """`self` should always be homogeneous."""
        for m in self.data:
            return self.deg_mon(m)

    # setters ----------------------------
    @classmethod
    def add_gen(cls, k: str, deg):
        """Add a new generator and return it."""
        cls._gen_names.append(k)
        cls._gen_degs.append(deg)
        m = (0,) * (len(cls._gen_names) - 1) + (1,)
        return cls(m).simplify()

    @classmethod
    def remove_gen(cls, k: str):
        """If `k`=0 in the algebra with relations simplified, call this function to remove this generator."""
        i = cls._gen_names.index(k)
        m_k = (0,) * i + (1,)
        del cls._rels[m_k]

        def f(_m):
            return _m[:i] + _m[i+1:]
        cls._rels = {f(_m): {f(m1) for m1 in _v} for _m, _v in cls._rels.items()}
        cls._gen_names = f(cls._gen_names)
        cls._gen_degs = f(cls._gen_degs)

    @classmethod
    def rename(cls, old_name, new_name):
        """Rename a generator."""
        i = cls._gen_names.index(old_name)
        cls._gen_names[i] = new_name

    @classmethod
    def add_gens(cls, names_degs):
        """Add generators. names_degs is a list of tuples (name, deg)."""
        for nd in names_degs:
            cls._gen_names.append(nd[0])
            cls._gen_degs.append(nd[1])

    @classmethod
    def add_rel(cls, rel):
        """Add a relation."""
        if not rel:
            return
        if type(rel) is not set:
            if not rel.is_homo():
                raise ValueError(f'relation {rel} not homogeneous!')
            hq = [(rel.deg(), rel.data)]
        else:
            if not cls(rel).is_homo():
                raise ValueError(f'relation {cls(rel)} not homogeneous!')
            hq = [(cls(rel).deg(), rel)]
        while hq:
            deg, r = heapq.heappop(hq)
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
                cls._rels[m] = r - {m}

    @classmethod
    def simplify_rels(cls):
        """Simplify `cls._rels`."""
        for m in cls._rels:
            cls._rels[m] = cls.simplify_data(cls._rels[m])

    @classmethod
    def simplify_data(cls, data: set):
        """Simplify `data` by relations."""
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
            cls.add_rel({x2, y})
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
    def get_lead(cls, data):
        return max(data, key=cls._key) if cls._key else max(data)

    @classmethod
    def gen(cls, k: str):
        """Return a generator."""
        i = cls._gen_names.index(k)
        m = (0,) * i + (1,)
        return cls(m).simplify() if cls.auto_simplify else cls(m)

    @classmethod
    def get_generators(cls, ideal=None):
        """Return a minimal generating set of `cls._rels` or `ideal`."""
        rel_gens = []
        if ideal is None:
            rels_backup = cls._rels
            cls._rels = {}
            try:
                for mon in sorted(rels_backup, key=cls.deg_mon):
                    rel_data = rels_backup[mon] | {mon}
                    if not cls.is_reducible(mon):
                        rel_gens.append(rel_data)
                        cls.add_rel(rel_data)
            finally:
                cls._rels = rels_backup
        else:
            rels_backup = cls._rels.copy()
            try:
                def deg_data(data):
                    for m in data:
                        return cls(m).deg()
                for rel_data in sorted(ideal, key=deg_data):
                    rel_data = cls.simplify_data(rel_data)
                    if rel_data:
                        rel_gens.append(rel_data)
                        cls.add_rel(rel_data)
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

    @classmethod
    def ann(cls, x):
        """Return the groebner basis for the ideal {y | xy=0}."""
        rels_backup = cls._rels.copy()
        key_backup = cls._key
        d = x.deg()
        try:
            X = cls.add_gen('X_{ann}', d)
            num_gen = len(X._gen_names)
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
            cls._key = key_backup
        result = []
        for m in rels:
            if len(m) == num_gen:
                result_i = cls.zero()
                for m1 in itertools.chain((m,), rels[m]):
                    result_i += cls(mymath.rstrip_tuple(m1[:-1])) * (x ** (m1[-1] - 1))
                result.append(result_i.data)
        return result

    def evaluation(self, image_gens):
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
    def print_tex(cls, show_gb=False):
        """For pdflatex."""
        print(f"Generators: ${', '.join(cls._gen_names)}$.\\\\")
        print(f"Degrees: ${', '.join(map(str, cls._gen_degs))}$")
        print("Relations:\\\\")
        if show_gb:
            for m in cls._rels:
                print(f"${cls(m)} = {cls(cls._rels[m])}$\\\\")
        else:
            for data in cls.get_generators():
                lead = cls.get_lead(data)
                print(f"${cls(lead)} = {cls(data - {lead})}$\\\\")

    @classmethod
    def display_alg(cls, show_gb=False):
        """For Jupyter notebook."""
        from IPython.display import HTML
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
            td = "\\begin{aligned}"
            for m in cls._rels:
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\"
            td += "\\end{aligned}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{aligned}"
            for data in cls.get_generators():
                lead = cls.get_lead(data)
                td += f"{cls(lead)} &= {cls(data - {lead})}\\\\"
            td += "\\end{aligned}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        result = "<table>" + tr1 + tr2 + "</table>" + "<table>" + tr3 + "</table>"
        return HTML(result)

    @classmethod
    def display_ideal(cls, gb: List[set]):
        from IPython.display import HTML
        result = "\\begin{aligned}"
        for data in gb:
            result += f"&{cls(data)}\\\\"
        result += "\\end{aligned}"
        return HTML(result)


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
               '_rels': {}, '_key': key, 'auto_simplify': True}
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
        from IPython.display import HTML
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
            td = "\\begin{aligned}"
            for m in cls._rels:
                td += f"{cls(m)} &= {cls(cls._rels[m])}\\\\"
            td += "\\end{aligned}"
        else:
            tr3 = "<th>Relations</th>"
            td = "\\begin{aligned}"
            for data in cls.get_generators():
                lead = cls.get_lead(data)
                td += f"{cls(lead)} &= {cls(data - {lead})}\\\\"
            td += "\\end{aligned}"
        td = f'<td>{td}</td>'
        tr3 += td
        tr3 = f"<tr>{tr3}</tr>"

        tr4 = '<th>Differentials</th>'
        td = '\\begin{aligned}'
        for g, dg in zip(cls._gen_names, cls._gen_diff):
            td += f'd({g}) &= {cls(dg)}\\\\'
        td += '\\end{aligned}'
        td = f'<td>{td}</td>'
        tr4 += td
        tr4 = f'<tr>{tr4}</tr>\n'

        result = '<table>\n' + tr1 + tr2 + '</table>' +\
                 '<table>' + tr3 + tr4 + '</table>'
        return HTML(result)


class SubRing:
    def __init__(self, d_max):
        self.d_max = d_max
        self.data = [[] for _ in range(d_max + 1)]

    def generate_non_comm(self, gens):
        if len(gens) == 0:
            raise ValueError("Need at least one generator")
        else:
            type_gen = type(gens[0])
            self.data[0].append([type_gen.unit(), type_gen.unit().get_mon()])
            for gen in gens:
                self.add_gen_non_comm(gen)

    def generate_comm(self, gens):
        if len(gens) == 0:
            raise ValueError("Need at least one generator")
        else:
            type_gen = type(gens[0])
            self.data[0].append([type_gen.unit(), type_gen.unit().get_mon()])
            for gen in gens:
                self.add_gen_comm(gen)

    def __str__(self):
        result = ""
        for d in range(self.d_max + 1):
            result += str(d) + ":\n"
            for base, _ in self.data[d]:
                result += str(base) + "\n"
        return result

    def add_gen_non_comm(self, gen):
        deg = gen.deg()
        if deg == 0:
            raise ValueError("the generator should have positive degree")
        for d in range(0, self.d_max - deg + 1):
            for base, _ in self.data[d]:
                new_base = base * gen
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

                new_base = gen * base  # non-commutative
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

    def add_gen_comm(self, gen):
        deg = gen.deg()
        if deg == 0:
            raise ValueError("the generator should have positive degree")
        for d in range(0, self.d_max - deg + 1):
            for base, _ in self.data[d]:
                new_base = base * gen
                for p, m in self.data[d + deg]:
                    if m in new_base.data:
                        new_base -= p
                if new_base:
                    self.data[d + deg].append([new_base, new_base.get_mon()])

    def basis(self, deg):
        return (gen for gen, _ in self.data[deg])


# noinspection PyUnresolvedReferences,PyArgumentList
class QuoRing:
    """ This class must be inherited with another ring """
    d_max = None
    ideal = None  # type: linalg.GradedVectorSpaceMod2

    def __mul__(self, other):
        prod = super().__mul__(other)
        res = self.ideal.res(prod)
        return res

    @classmethod
    def init(cls, d_max):
        cls.d_max = d_max
        cls.ideal = linalg.GradedVectorSpaceMod2(d_max)

    @classmethod
    def add_relations(cls, rels, get_mon=None):
        for rel in rels:
            if rel:
                deg = rel.deg()
                for d in range(cls.d_max - deg + 1):
                    cls.ideal.add_vectors((super(QuoRing, a).__mul__(rel) for a in super().basis(d)), deg + d, get_mon)

    @classmethod
    def basis_mons(cls, deg):
        for m in super().basis_mons(deg):
            if m not in cls.ideal.get_mons(deg):
                yield m

    @classmethod
    def basis(cls, deg):
        return (cls(m) for m in cls.basis_mons(deg))


class FreeModule:
    type_ring = None

    def __init__(self, data: dict):
        if type(data) is dict:
            self.data = dict((mon, coeff) for mon, coeff in data.items() if coeff)
        else:
            raise TypeError

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        data = copy.deepcopy(self.data)
        for key in other.data:
            if key in data:
                data[key] += other.data[key]
            else:
                data[key] = other.data[key].copy()
        return type(self)(data)

    def __iadd__(self, other):
        data = self.data
        for key in other.data:
            if key in data:
                data[key] += other.data[key]
            else:
                data[key] = other.data[key].copy()
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __sub__(self, other):
        data = copy.deepcopy(self.data)
        for key in other.data:
            if key in data:
                data[key] -= other.data[key]
            else:
                data[key] = -other.data[key]
        return type(self)(data)

    def __isub__(self, other):
        data = self.data
        for key in other.data:
            if key in data:
                data[key] -= other.data[key]
            else:
                data[key] = -other.data[key]
        self.data = dict((mon, coeff) for mon, coeff in self.data.items() if coeff)
        return self

    def __rmul__(self, other):
        if type(other) is self.type_ring:
            data = {}
            for key in self.data:
                data[key] = other * self.data[key]
            return type(self)(data)
        else:
            return NotImplemented

    def __str__(self):
        result = ""
        for key in sorted(self.data):
            if result != "":
                result += "+"
            str_coeff = str(self.data[key])
            if '+' in str_coeff:
                result += "({}){}".format(str_coeff, key)
            elif str_coeff == "1":
                result += str(key)
            else:
                result += "{}{}".format(str_coeff, key)
        if result == "":
            result += "0"
        return result

    @classmethod
    def set_ring(cls, type_ring):
        cls.type_ring = type_ring

    @classmethod
    def gen(cls, key):
        """ return a generator of type FM """
        return cls({key: cls.type_ring.unit()})

    def is_zero(self) -> bool:
        return len(self.data) == 0

    def get_mon(self) -> tuple:
        key = max(self.data)
        return key, self.data[key].get_mon()

    def has_mon(self, mon: tuple) -> bool:
        if mon[0] in self.data:
            return mon[1] in self.data[mon[0]].data
        else:
            return False


class FreeModuleMod2:
    type_ring = None

    def __init__(self, data: Union[set, Tuple[tuple, Union[str, int]]]):
        if type(data) is set:
            self.data = data
        elif type(data) is tuple:  # monomial
            self.data = {data}  # type: Set[Tuple[tuple, Union[str, int]]]
        else:
            raise TypeError("{} can not initialize {}".format(data, type(self).__name__))

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        return type(self)(self.data ^ other.data)

    def __iadd__(self, other):
        self.data ^= other.data
        return self

    def __sub__(self, other):
        return type(self)(self.data ^ other.data)

    def __isub__(self, other):
        self.data ^= other.data
        return self

    def __rmul__(self, other):
        if type(other) is self.type_ring:
            data = set()
            for mon in self.data:
                prod = other * self.type_ring(mon[0])
                data.symmetric_difference_update((m, mon[1]) for m in prod.data)
            return type(self)(data)
        else:
            return NotImplemented

    def __str__(self):
        result = ""
        for m in self.data:
            if result != "":
                result += "+"
            str_coeff = str(self.type_ring(m[0]))
            if str_coeff == "1":
                result += str(m[1])
            else:
                result += "{}{}".format(str_coeff, m[1])
        if result == "":
            result += "0"
        return result

    @classmethod
    def set_ring(cls, type_ring):
        cls.type_ring = type_ring

    @classmethod
    def gen(cls, key):
        """ return a generator of type FM """
        return cls(((), key))

    @classmethod
    def zero(cls):
        return cls(set())

    def is_zero(self) -> bool:
        return len(self.data) == 0

    def get_mon(self) -> tuple:
        return max(self.data)

    def has_mon(self, mon: tuple) -> bool:
        return mon in self.data


def alg_may(n_max):
    R = GbAlgMod2.new_alg()
    gens = []

    def b(_i, _j):
        return R.gen(f"b^{_i}_{_j}") if _j > 1 else R.gen(f"h_{_i}") ** 2

    # add generators
    for s in range(n_max):
        gens.append((f"h_{s}", 1, (1, -2 ** s)))
    for t in range(1, n_max + 1):
        for s in reversed(range(0, t)):
            if t - s > 1:
                gens.append((f"b^{s}_{t-s}", 2, (0, -2 * (2 ** t - 2 ** s))))
    for s in range(n_max - 2):
        gens.append((f"h_{s}(1)", 2, (2, -9 * 2 ** s)))
    for s in range(n_max - 4):
        gens.append((f"h_{s}(1,2)", 3, (3, -49 * 2 ** s)))
    for s in range(n_max - 4):
        gens.append((f"h_{s}(1,3)", 3, (3, -41 * 2 ** s)))
    gens.sort(key=lambda _x: _x[2], reverse=True)
    R.add_gens(gens)

    # add relations
    for j in range(n_max + 1):
        for s in reversed(range(0, n_max - j + 1)):
            t = s + j
            rel = sum((b(s, k-s) * b(k, t-k) for k in range(s + 1, t)), R.zero())
            R.add_rel(rel)
            # print(s, t)
    for s in range(n_max - 2):
        rel = R.gen(f"h_{s}(1)") * R.gen(f"h_{s}(1)") + b(s, 3) * b(s+1, 1) + \
              b(s, 2) * b(s+1, 2)
        R.add_rel(rel)
        # print(s)
    R.reduce_frob()
    R.print_tex()
    return R


# 140, 248, 283, 415, 436, 612, 600, 588, 701, 782
