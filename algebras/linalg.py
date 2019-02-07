"""A module for linear algebra.

This module provides interfaces for algebras.
Warning: If you pass a set as a parameter to a method in this module,
it might be modified. If you pass an algebra instead, the method would create
a shallow copy of it.
"""

from typing import Iterable, Set, Tuple, List, Callable, Union
import copy


_t_mon = Union[tuple, frozenset, int, str]  # todo: generic type


class VectorSpaceMod2:
    """This class is for modeling vector spaces."""
    def __init__(self, vectors=None, *, data=None, key=None):
        self.key = key
        self.data = data if data is not None else []
        if vectors is not None:
            self.add_vectors(vectors)

    def copy(self):
        """Return a deep copy."""
        return VectorSpaceMod2(data=copy.deepcopy(self.data))

    def add_vectors(self, vectors: Iterable):
        """Add vectors to this vector space."""
        for v in vectors:
            self.add_v(v)

    def add_v(self, v):
        if type(v) is not set:
            v = v.data.copy()
        for w, mw in self.data:
            if mw in v:
                v ^= w
        if v:
            self.data.append((v, max(v, key=self.key) if self.key else max(v)))

    def simplify(self) -> "VectorSpaceMod2":
        for i in range(len(self.data) - 1, 0, -1):
            v, mv = self.data[i]
            for j in range(i):
                w = self.data[j][0]
                if mv in w:
                    w ^= v
        return self

    def res(self, vector):
        """Return vector mod self."""
        v = vector if type(vector) is set else vector.data.copy()
        for w, mw in self.data:
            if mw in v:
                v ^= w
        return v if type(vector) is set else type(vector)(v)

    def quotient(self, other: "VectorSpaceMod2") -> "VectorSpaceMod2":
        """Return self/other assuming other is a subspace."""
        result = other.data.copy()
        n = len(result)
        for v, mv in self.data:
            v = v.copy()
            for v1, mv1 in result:
                if mv1 in v:
                    v ^= v1
            if v:
                result.append((v, max(v, key=self.key) if self.key else max(v)))
        return VectorSpaceMod2(data=result[n:])

    def get_mons(self) -> Set[_t_mon]:
        """Return the leading monomials."""
        return set(m for _, m in self.data)

    def basis(self, type_alg) -> Iterable:
        return (type_alg(v) for v, mv in self.data)

    def dim(self) -> int:
        """Return the dimension of the vector space."""
        return len(self.data)


class GradedVectorSpaceMod2:
    """ a graded version of VectorSpaceMod2 """
    def __init__(self, d_max, *, key=None):
        self.d_max = d_max
        self.data = [VectorSpaceMod2(key=key) for _ in range(d_max + 1)]  # type: List[VectorSpaceMod2]

    def add_vectors(self, vectors: Iterable[_t_mon], deg: int):
        """ add vectors to this vector space """
        assert(0 <= deg <= self.d_max)
        self.data[deg].add_vectors(vectors)

    def res(self, vector, is_homogeneous=True):
        """ return vector mod this VectorSpace"""
        if not vector:
            return vector
        assert(vector.deg() <= self.d_max)
        if is_homogeneous:
            return self.data[vector.deg()].res(vector)
        else:
            degs = set(vector.deg_mon(mon) for mon in vector.data)
            res = vector
            for deg in degs:
                res = self.data[deg].res(vector)
            return res

    def get_mons(self, deg):
        """ return the leading monomials """
        return self.data[deg].get_mons()


class LinearMapMod2:  # TODO: add key
    """ this class is for modeling linear maps f: V leftrightarrow W: g"""
    def __init__(self):
        # self.maps is for f:V->W and self.inv_maps is for g:W->V
        self.maps = []  # type: List[Tuple[Set[_t_mon], _t_mon, Set[_t_mon]]]
        self.inv_maps = []  # type: List[Tuple[Set[_t_mon], _t_mon, Set[_t_mon]]]
        self.kernel = []  # type: List[Tuple[Set[_t_mon], _t_mon]]

    def clear(self):
        self.__init__()

    def add_maps(self, maps: Iterable):
        for src, tgt in maps:
            v = src.data.copy()  # type: Set[_t_mon]
            fv = tgt.data.copy()  # type: Set[_t_mon]
            for v1, mv1, fv1 in self.maps:
                if mv1 in v:
                    v ^= v1
                    fv ^= fv1
            if not v:
                if fv:
                    raise ValueError("linear map not well-defined")
            else:
                self.maps.append((v, max(v), fv))

                w = fv.copy()
                gw = v.copy()
                for w1, mw1, gw1 in self.inv_maps:
                    if mw1 in w:
                        w ^= w1
                        gw ^= gw1
                if not w:  # we get gw in the kernel
                    for v1, mv1 in self.kernel:
                        if mv1 in gw:
                            gw ^= v1
                    self.kernel.append((gw, max(gw)))
                else:
                    self.inv_maps.append((w, max(w), gw))

    def f(self, vector):
        """ return f(vector) """
        type_vector = type(vector)
        v = vector if type_vector is set else vector.data.copy()
        result = set()
        for v1, mv1, fv1 in self.maps:
            if mv1 in v:
                v ^= v1
                result ^= fv1
        if not v:
            return result if type_vector is set else type_vector(result)
        else:
            return None

    def g(self, vector):
        """ return f^{-1}(vector) """
        w = vector if type(vector) is set else vector.data.copy()
        result = set()
        for w1, mw1, gw1 in self.inv_maps:
            if mw1 in w:
                w ^= w1
                result ^= gw1
        if not w:
            return type(vector)(result)
        else:
            return None

    def present_kernel(self, type_vector):
        for v, mv in self.kernel:
            print(type_vector(v))


class LinearMapKernelMod2:
    """ this is a optimized version of LinearMapMod2 that focus on computing the kernel """
    def __init__(self, *, key=None):
        # self.inv_maps is for g:W->V
        self.key = key
        self.inv_maps = []  # type: List[Tuple[Set[_t_mon], _t_mon, Set[_t_mon]]]
        self.kernel = VectorSpaceMod2(key=key)

    def clear(self):
        self.__init__()

    def add_maps(self, maps: Iterable[tuple]):
        for src, tgt in maps:
            gw = src if type(src) is set else src.data.copy()  # type: Set[_t_mon]
            w = tgt if type(tgt) is set else tgt.data.copy()  # type: Set[_t_mon]
            for w1, mw1, gw1 in self.inv_maps:
                if mw1 in w:
                    w ^= w1
                    gw ^= gw1
            if not w:  # we get gw in the kernel
                self.kernel.add_v(gw)
            else:
                self.inv_maps.append((w, max(w, key=self.key) if self.key else max(w), gw))

    def g(self, vector):
        """ return f^{-1}(vector) """
        w = vector if type(vector) is set else vector.data.copy()
        result = set()
        for w1, mw1, gw1 in self.inv_maps:
            if mw1 in w:
                w ^= w1
                result ^= gw1
        if not w:
            return type(vector)(result)
        else:
            return None

    def image(self) -> VectorSpaceMod2:
        """ warning: the return should not be modified """
        return VectorSpaceMod2(data=[(w, mw) for w, mw, gw in self.inv_maps], key=self.key)
