import pickle
from notations import *

x = Sq(4) * Sq(4) * Sq(4)
print(x)

f = open("tmp.pickle", "wb")
pickle.dump(x, f)
f.close()

f = open("tmp.pickle", "rb")
y = pickle.load(f)
f.close()

print(y)
