from algebras.operations import DualSteenrod
import itertools


class CobarSteenrod:
    def __init__(self, data: set):
        self.data = data

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        result = ""
        for m in self.data:
            if len(result) > 0:
                result += "+"
            result += self.str_mon(m)
        if result == "":
            result += "0"
        return result

    def __add__(self, other):
        return type(self)(self.data ^ other.data)

    def __mul__(self, other):
        data = set()
        for m1 in self.data:
            for m2 in other.data:
                data ^= {m1 + m2}
        return type(self)(data)

    @classmethod
    def tensor(cls, *args: DualSteenrod):
        multi_data = tuple(r.data for r in args)
        iter_prod = itertools.product(*multi_data)
        data = set(m for m in iter_prod)
        return cls(data)

    @staticmethod
    def str_mon(mon):
        if len(mon) == 0:
            return "1"
        result = "["
        for i in range(len(mon)):
            result += DualSteenrod.str_mon(mon[i])
            if i < len(mon) - 1:
                result += "|"
        result += "]"
        return result

    @classmethod
    def zero(cls):
        return cls(set())

    @classmethod
    def unit(cls):
        return cls({()})

    def diff(self):
        data = set()
        for mon in self.data:
            for i in range(len(mon)):
                m = mon[i]
                coprod = DualSteenrod(m).coprod()
                for m1, m2 in coprod.data:
                    if m1 and m2:
                        data ^= {mon[:i] + (m1, m2) + mon[i+1:]}
        return type(self)(data)
