import time
from collections import Counter
import array

N = 60_000_000
t0 = time.time()

# Memoized total stopping time using a cache array for the whole range
# reached during trajectories (grows beyond N, so use a dict for the
# overflow above a fixed array cache).
CACHE_SIZE = N * 4
cache = array.array('i', [0]) * 0  # placeholder, will build properly
cache = array.array('l', [-1]) * CACHE_SIZE  # -1 = unknown, sized array

def steps(n, cache=cache, CACHE_SIZE=CACHE_SIZE):
    path = []
    x = n
    while True:
        if x == 1:
            base = 0
            break
        if x < CACHE_SIZE and cache[x] != -1:
            base = cache[x]
            break
        path.append(x)
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
    # unwind, assigning steps counts
    total = base
    for v in reversed(path):
        total += 1
        if v < CACHE_SIZE:
            cache[v] = total
    return total

def early_merge_time(n, s):
    x, y = n, n + 1
    for t in range(1, s):
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        y = y >> 1 if y & 1 == 0 else 3 * y + 1
        if x == y:
            return t
    return None

agree = 0
early_merge = 0
merge_times = []
exceptions = []

for n in range(1, N + 1):
    sX = steps(n)
    sY = steps(n + 1)
    if sX == sY:
        agree += 1
        mt = early_merge_time(n, sX)
        if mt is not None:
            early_merge += 1
            merge_times.append(mt)
        else:
            exceptions.append(n)

print("total time:", time.time() - t0)
print(f"\nN={N}")
print(f"raw agreement rate: {agree/N:.5f}  ({agree}/{N})")
print(f"early-merge fraction: {early_merge/agree:.6f}  ({early_merge}/{agree})")
print(f"EXCEPTIONS (agree but no early merge): {len(exceptions)}")
if exceptions:
    print("first few exceptions:", exceptions[:20])

if merge_times:
    print(f"\nmean early-merge time: {sum(merge_times)/len(merge_times):.4f}")
    print(f"max early-merge time: {max(merge_times)}")

with open('/home/claude/exceptions.txt', 'w') as f:
    for e in exceptions:
        f.write(f"{e}\n")
