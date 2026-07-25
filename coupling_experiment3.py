import time
from collections import Counter

def total_stopping_time(n):
    """steps(n): number of 3x+1/x/2 steps until reaching 1 (stops at 1)."""
    t = 0
    x = n
    while x != 1:
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        t += 1
    return t

def first_merge_time(n, cap):
    """
    tau_couple: first t such that T^t(n) == T^t(n+1), applying the RAW map
    (no special-casing at 1), but we cap the search at 'cap' steps
    (cap = min(steps(n), steps(n+1)) + small buffer) so we don't pick up
    spurious coincidences from both trajectories cycling in 1->4->2->1
    after having already reached 1 independently.
    Returns t if merged within cap, else None.
    """
    x, y = n, n + 1
    for t in range(1, cap + 1):
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        y = y >> 1 if y & 1 == 0 else 3 * y + 1
        if x == y:
            return t
    return None

N = 300_000
t0 = time.time()

agree = 0          # steps(n) == steps(n+1)
literal_merge = 0  # AND they literally coincide in value before/at reaching 1
merge_times = []
steps_cache_prev = None

# compute steps(n) for n=1..N+1 once, reuse
steps_arr = [0] * (N + 2)
for n in range(1, N + 2):
    steps_arr[n] = total_stopping_time(n)

for n in range(1, N + 1):
    sX, sY = steps_arr[n], steps_arr[n + 1]
    if sX == sY:
        agree += 1
        cap = sX + 2  # steps to reach 1 (should merge by then if ever)
        mt = first_merge_time(n, cap)
        if mt is not None:
            literal_merge += 1
            merge_times.append(mt)

print("time:", time.time() - t0)
print(f"N={N}")
print(f"raw agreement rate P(steps(n)=steps(n+1)): {agree/N:.5f}  ({agree}/{N})")
print(f"of those, literal-merge fraction: {literal_merge/agree:.5f}  ({literal_merge}/{agree})")
print(f"literal-merge as fraction of ALL n: {literal_merge/N:.5f}")

print(f"\nmean coupling time (given merge): {sum(merge_times)/len(merge_times):.4f}")
print(f"max coupling time: {max(merge_times)}")
print(f"min coupling time: {min(merge_times)}")

c = Counter(merge_times)
max_t = max(merge_times)
total_merged = len(merge_times)
cum = 0
survival = []
for t in range(0, max_t + 2):
    cum += c.get(t, 0)
    survival.append((total_merged - cum) / total_merged)

print("\nt : P(tau_couple > t | merged)")
for t in [0,1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,70,100]:
    if t < len(survival):
        print(f"{t:4d} : {survival[t]:.6f}")

with open('/home/claude/merge_times.txt', 'w') as f:
    for t in merge_times:
        f.write(f"{t}\n")
