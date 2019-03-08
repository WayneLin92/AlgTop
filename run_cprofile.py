import cProfile
import pstats
from pstats import SortKey
from algebras.constructions import alg_bij

pr = cProfile.Profile()
pr.enable()

# run
alg_bij(7)

pr.disable()
with open("log_cProfile_alg_bij.txt", "w") as file:
    sort_by = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=file).strip_dirs().sort_stats(sort_by)
    ps.print_stats()
