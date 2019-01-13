import cProfile
import io
import pstats
from pstats import SortKey
from spec.resolution import test_ext

pr = cProfile.Profile()
pr.enable()

# run
test_ext()

pr.disable()
with open("log_cProfile.txt", "w") as file:
    sort_by = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=file).strip_dirs().sort_stats(sort_by)
    ps.print_stats()
