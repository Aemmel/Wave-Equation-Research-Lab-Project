"""
What to test?
test multiple versions of the program:
    - no parallelization (standard way of doing it)
    - parallelized with fine grains (write "#pragma omp parallel for" very often)
    - coarse grained (only one "#pragma omp parallel" block)
  -> has to be recompiled. Can't be called from here
what parameters to test:
    - cores (1...8)
    - grid size (100x100 to max my computer can handle) for different core number
    - chunk size for parallel loops (fixed grid)
    - 
"""

from subprocess import run
from time import perf_counter



grid_sizes = list(range(500, 4200, 200)) # 4000x4000 seems to still work. 5000x5000 crashes
cores = list(range(1, 9))


begin = perf_counter()
run("./build/wave")
end = perf_counter()

print("took: " + str(end-begin))