"""A module for linear algebra.

This module provides interfaces for algebras.
Warning: If you pass a set as a parameter to a method in this module,
it might be modified. If you pass an algebra instead, the method would create
a shallow copy of it.
"""
import copy
import operator
import collections
from typing import TypeVar, Tuple, List, Set, Dict, Iterable, Any
from algebras import myerror

_t_mon = TypeVar('_t_mon')
_t_v = Set[_t_mon]
_t_data = List[Tuple[_t_v, _t_mon]]


class VectorSpaceMod2:
    """A class for vector spaces."""
    def __init__(self, vectors=None, *, data=None, key=None):
        self.key = key
        self.data = data or []
        if vectors is not None:
            self.add_vectors(vectors)

    def copy(self):
        """Return a deep copy."""
        return VectorSpaceMod2(data=copy.deepcopy(self.data))

    # setters ----------------
    def add_v_set(self, v):
        """Add a single vector efficiently."""
        for v1, mv1 in self.data:
            if mv1 in v:
                v ^= v1
        if v:
            self.data.append((v, max(v, key=self.key) if self.key else max(v)))

    def add_vectors_set(self, vectors: Iterable[_t_v]):
        """Add vectors efficiently."""
        for v in vectors:
            for v1, mv1 in self.data:
                if mv1 in v:
                    v ^= v1
            if v:
                self.data.append((v, max(v, key=self.key) if self.key else max(v)))

    def add_v(self, v):
        """Add a single vector."""
        if type(v) is not set:
            v = v.data.copy()
        self.add_v_set(v)

    def add_vectors(self, vectors: Iterable):
        """Add vectors."""
        for v in vectors:
            self.add_v(v)

    # modifiers ------------------
    def simplify(self) -> "VectorSpaceMod2":
        """Simplify the basis such that it forms a block matrix (I, A)."""
        for i in range(len(self.data) - 1, 0, -1):
            v, mv = self.data[i]
            for j in range(i):
                w = self.data[j][0]
                if mv in w:
                    w ^= v
        return self

    # getters ------------------
    def get_mons(self):
        """Return the leading monomials."""
        return map(operator.itemgetter(1), self.data)

    def basis(self, type_alg=set):
        """Return a basis of the vector space."""
        vectors = map(operator.itemgetter(0), self.data)
        return vectors if type_alg is set else map(type_alg, vectors)

    def present(self, type_alg=set):
        print("Vector Space:")
        for r in self.basis(type_alg):
            print(r)

    @property
    def dim(self) -> int:
        """Return the dimension of the vector space."""
        return len(self.data)

    # functions -----------------
    def res_set(self, v: set):
        """Return v mod self efficiently."""
        for v1, mv1 in self.data:
            if mv1 in v:
                v ^= v1
        return v

    def res(self, vector):
        """Return vector mod self."""
        v = vector if type(vector) is set else vector.data.copy()
        self.res_set(v)
        return v if type(vector) is set else type(vector)(v)

    def __truediv__(self, other: "VectorSpaceMod2"):
        """Return the quotient space self/other."""
        result = VectorSpaceMod2()
        result.add_vectors_set(map(other.res_set, self.basis(set)))
        return result

    def __iadd__(self, other: "VectorSpaceMod2"):
        """Expand self by another vector space."""
        self.add_vectors_set(other.basis(set))


class GradedVectorSpaceMod2:
    """A graded version of VectorSpaceMod2."""
    def __init__(self, *, key=None):
        self.key = key
        self.data = {}  # type: Dict[Any, _t_data]

    # setters ----------------
    def add_v_set(self, v: _t_v, deg):
        """Add a single vector efficiently."""
        if deg in self.data:
            for v1, mv1 in self.data[deg]:
                if mv1 in v:
                    v ^= v1
            if v:
                self.data[deg].append((v, max(v, key=self.key) if self.key else max(v)))
        else:
            self.data[deg] = [(v, max(v, key=self.key) if self.key else max(v))]

    def add_vectors_set(self, vectors: Iterable[_t_v], deg):
        """Add vectors efficiently."""
        if deg not in self.data:
            self.data[deg] = []
        for v in vectors:
            for v1, mv1 in self.data[deg]:
                if mv1 in v:
                    v ^= v1
            if v:
                self.data[deg].append((v, max(v, key=self.key) if self.key else max(v)))

    def add_v(self, v):
        """Add a single vector from an algebra."""
        deg = v.deg()
        v = v.data.copy()
        self.add_v_set(v, deg)

    def add_vectors(self, vectors):
        """Add vectors from an algebra."""
        for v in vectors:
            self.add_v(v)

    # modifiers ------------------
    def simplify(self) -> "GradedVectorSpaceMod2":
        """Simplify the basis such that it forms a block matrix (I, A)."""
        for deg in self.data:
            for i in range(len(self.data[deg]) - 1, 0, -1):
                v, mv = self.data[deg][i]
                for j in range(i):
                    w = self.data[deg][j][0]
                    if mv in w:
                        w ^= v
        return self

    # getters ------------------
    def get_mons(self, deg):
        """Return the leading monomials."""
        return map(operator.itemgetter(1), self.data[deg]) if deg in self.data else ()

    def basis(self, deg, type_alg):
        """Return a basis of the vector space."""
        return map(type_alg, map(operator.itemgetter(0), self.data[deg])) if deg in self.data else ()

    def dim(self, deg) -> int:
        """Return the dimension of the vector space."""
        return len(self.data[deg]) if deg in self.data else 0

    # functions -----------------
    def res_set(self, v, deg):
        """Return v mod self efficiently."""
        if deg in self.data:
            for v1, mv1 in self.data[deg]:
                if mv1 in v:
                    v ^= v1
        return v

    def res(self, vector):
        """ return vector mod this VectorSpace"""
        deg = vector.deg()
        v = vector.data.copy()
        self.res_set(v, deg)
        return type(vector)(v)


class LinearMapMod2:
    """Linear map f: V leftrightarrow W: g."""
    def __init__(self, *, key=None):
        self.key = key
        self._domain = []  # type: _t_data
        self._f = []  # type: List[set]
        self._image = []  # type: _t_data
        self._g = []  # type: List[set]
        self._kernel = VectorSpaceMod2(key=key)

    # setters ----------------
    def add_maps_set(self, maps: Iterable[Tuple[set, set]]):
        """Add maps efficiently."""
        for v, fv in maps:
            for vm1, fv1 in zip(self._domain, self._f):
                if vm1[1] in v:
                    v ^= vm1[0]
                    fv ^= fv1
            if not v:
                if fv:
                    raise myerror.MyValueError("Incompatible linear map.")
            else:
                self._domain.append((v, max(v, key=self.key) if self.key else max(v)))
                self._f.append(fv)

                w = fv.copy()
                gw = v.copy()
                for wm1, gw1 in zip(self._image, self._g):
                    if wm1[1] in w:
                        w ^= wm1[0]
                        gw ^= gw1
                if not w:
                    self._kernel.add_v_set(gw)
                else:
                    self._image.append((w, max(w, key=self.key) if self.key else max(w)))
                    self._g.append(gw)

    def add_maps(self, maps):
        """Add maps."""
        self.add_maps_set((v.data.copy(), fv.data.copy()) for v, fv in maps)

    def add_map(self, v, fv):
        self.add_maps(((v, fv),))

    # getters ------------------
    @property
    def domain(self):
        return VectorSpaceMod2(data=self._domain)

    @property
    def image(self):
        return VectorSpaceMod2(data=self._image)

    @property
    def kernel(self):
        return self._kernel

    # functions -----------------
    def f(self, vector):
        """Return f(vector)."""
        v = vector.data.copy()
        result = set()
        for vm1, fv1 in zip(self._domain, self._f):
            if vm1[1] in v:
                v ^= vm1[0]
                result ^= fv1
        return None if v else type(vector)(result)

    def g(self, vector):
        """Return f^{-1}(vector)."""
        w = vector.data.copy()
        result = set()
        for wm1, gw1 in zip(self._image, self._g):
            if wm1[1] in w:
                w ^= wm1[0]
                result ^= gw1
        return None if w else type(vector)(result)


class LinearMapKMod2:
    """This is an optimized version of LinearMapMod2 that focus on computing the kernel.

    Warning: Incompatible maps will not be reported."""
    def __init__(self, *, key=None):
        self.key = key
        self._image = []  # type: _t_data
        self._g = []  # type: List[set]
        self._kernel = VectorSpaceMod2(key=key)

    # setters ----------------
    def add_maps_set(self, maps: Iterable[Tuple[set, set]]):
        """Add maps efficiently."""
        for gw, w in maps:
            for wm1, gw1 in zip(self._image, self._g):
                if wm1[1] in w:
                    w ^= wm1[0]
                    gw ^= gw1
            if not w:
                self._kernel.add_v_set(gw)
            else:
                self._image.append((w, max(w, key=self.key) if self.key else max(w)))
                self._g.append(gw)

    def add_maps(self, maps):
        """Add maps."""
        self.add_maps_set((v.data.copy(), fv.data.copy()) for v, fv in maps)

    def add_map(self, v, fv):
        self.add_maps(((v, fv),))

    # getters ------------------
    @property
    def image(self):
        return VectorSpaceMod2(data=self._image)

    @property
    def kernel(self):
        return self._kernel

    def g(self, vector):
        """Return f^{-1}(vector)."""
        w = vector.data.copy()
        result = set()
        for wm1, gw1 in zip(self._image, self._g):
            if wm1[1] in w:
                w ^= wm1[0]
                result ^= gw1
        return None if w else type(vector)(result)


class GradedLinearMapMod2:
    """A graded version of LinearMapMod2."""
    MyTuple = collections.namedtuple('MyTuple', ('domain', 'f', 'image', 'g', 'kernel'))

    def __init__(self, *, key=None):
        self.key = key
        self.data = {}  # type: Dict[Any, GradedLinearMapMod2.MyTuple]

    # setters ----------------
    def add_maps_set(self, maps: Iterable[Tuple[set, set]], deg):
        """Add maps efficiently."""
        if deg not in self.data:
            self.data[deg] = self.MyTuple([], [], [], [], VectorSpaceMod2(key=self.key))
        linmap = self.data[deg]
        for v, fv in maps:
            for vm1, fv1 in zip(linmap.domain, linmap.f):
                if vm1[1] in v:
                    v ^= vm1[0]
                    fv ^= fv1
            if not v:
                if fv:
                    raise myerror.MyValueError("Incompatible linear map.")
            else:
                linmap.domain.append((v, max(v, key=self.key) if self.key else max(v)))
                linmap.f.append(fv)

                w = fv.copy()
                gw = v.copy()
                for wm1, gw1 in zip(linmap.image, linmap.g):
                    if wm1[1] in w:
                        w ^= wm1[0]
                        gw ^= gw1
                if not w:
                    linmap.kernel.add_v_set(gw)
                else:
                    linmap.image.append((w, max(w, key=self.key) if self.key else max(w)))
                    linmap.g.append(gw)

    def add_maps(self, maps, deg):
        """Add maps."""
        self.add_maps_set(((v.data.copy(), fv.data.copy()) for v, fv in maps), deg)

    def add_map(self, v, fv):
        deg = v.deg()
        self.add_maps(((v, fv),), deg)

    # getters ------------------
    def domain(self, deg):
        return VectorSpaceMod2(data=self.data[deg].domain if deg in self.data else None)

    def image(self, deg):
        return VectorSpaceMod2(data=self.data[deg].image if deg in self.data else None)

    def kernel(self, deg):
        return self.data[deg].kernel if deg in self.data else VectorSpaceMod2()

    # functions -----------------
    def f(self, vector):
        """Return f(vector)."""
        deg = vector.deg()
        if deg in self.data:
            linmap = self.data[deg]
        else:
            return None
        v = vector.data.copy()
        result = set()
        for vm1, fv1 in zip(linmap.domain, linmap.f):
            if vm1[1] in v:
                v ^= vm1[0]
                result ^= fv1
        return None if v else type(vector)(result)

    def g(self, vector):
        """Return f^{-1}(vector)."""
        deg = vector.deg()
        if deg in self.data:
            linmap = self.data[deg]
        else:
            return None
        w = vector.data.copy()
        result = set()
        for wm1, gw1 in zip(linmap.image, linmap.g):
            if wm1[1] in w:
                w ^= wm1[0]
                result ^= gw1
        return None if w else type(vector)(result)


class GradedLinearMapKMod2:
    """A graded version of LinearMapKMod2."""
    MyTupleK = collections.namedtuple('MyTupleK', ('image', 'g', 'kernel'))

    def __init__(self, *, key=None):
        self.key = key
        self.data = {}  # type: Dict[Any, GradedLinearMapKMod2.MyTupleK]

    # setters ----------------
    def add_maps_set(self, maps: Iterable[Tuple[set, set]], deg):
        """Add maps efficiently."""
        if deg not in self.data:
            self.data[deg] = self.MyTupleK([], [], VectorSpaceMod2(key=self.key))
        linmap = self.data[deg]
        for gw, w in maps:
            for wm1, gw1 in zip(linmap.image, linmap.g):
                if wm1[1] in w:
                    w ^= wm1[0]
                    gw ^= gw1
            if not w:
                linmap.kernel.add_v_set(gw)
            else:
                linmap.image.append((w, max(w, key=self.key) if self.key else max(w)))
                linmap.g.append(gw)

    def add_maps(self, maps, deg):
        """Add maps."""
        self.add_maps_set(((v.data.copy(), fv.data.copy()) for v, fv in maps), deg)

    def add_map(self, v, fv):
        deg = v.deg()
        self.add_maps(((v, fv),), deg)

    # getters ------------------
    def image(self, deg):
        return VectorSpaceMod2(data=self.data[deg].image if deg in self.data else None)

    def kernel(self, deg):
        return self.data[deg].kernel if deg in self.data else VectorSpaceMod2()

    # functions -----------------
    def g(self, vector):
        """Return f^{-1}(vector)."""
        deg = vector.deg()
        if deg in self.data:
            linmap = self.data[deg]
        else:
            return None
        w = vector.data.copy()
        result = set()
        for wm1, gw1 in zip(linmap.image, linmap.g):
            if wm1[1] in w:
                w ^= wm1[0]
                result ^= gw1
        return None if w else type(vector)(result)

# 226, 302, 311, 412, 406
