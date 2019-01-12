from notations import *
from algebras.my_dyerlashof import *


# noinspection PyPep8Naming
def sQ(*I):
    result = MyDyerLashof.unit()
    for i in I:
        result *= MyDyerLashof.gen(i)
    return result
