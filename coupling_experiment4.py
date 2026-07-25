import time
from collections import Counter

def total_stopping_time(n):
    t = 0
    x = n
    while x != 1:
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        t += 1
    return t

def early_merge_time(n, s):
    """
    Search for t in [1, s-1] (STRICTLY before both reach 1) such that
    T^t(n) == T^t(n+1). Returns t if found (genuine early merge),
    else None (they only coincide trivially at t=s, both equal to 1).
    """
    x, y = n, n + 1
    for t in range(1, s):
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        y = y >> 1 if y & 1 == 0 else 3 * y + 1
        if x == y:
            return t
    return None

N = 5_000_000
t0 = time.time()

steps_arr = [0] * (N + 2)
for n in range(1, N + 2):
    steps_arr[n] = total_stopping_time(n)
print("steps computed:", time.time() - t0)

agree = 0
early_merge = 0
merge_times = []

for n in range(1, N + 1):
    sX, sY = steps_arr[n], steps_arr[n + 1]
    if sX == sY:
        agree += 1
        mt = early_merge_time(n, sX)
        if mt is not None:
            early_merge += 1
            merge_times.append(mt)

print("total time:", time.time() - t0)
print(f"\nN={N}")
print(f"raw agreement rate P(steps(n)=steps(n+1)): {agree/N:.5f}  ({agree}/{N})")
print(f"of agreeing pairs, EARLY-merge fraction (merge strictly before reaching 1): {early_merge/agree:.5f}  ({early_merge}/{agree})")
print(f"pairs that agree but ONLY coincide trivially at final value 1: {(agree-early_merge)/agree:.5f}  ({agree-early_merge}/{agree})")

if merge_times:
    print(f"\nmean early-merge time: {sum(merge_times)/len(merge_times):.4f}")
    print(f"max early-merge time: {max(merge_times)}")
    print(f"min early-merge time: {min(merge_times)}")

    c = Counter(merge_times)
    max_t = max(merge_times)
    total_merged = len(merge_times)
    cum = 0
    survival = []
    for t in range(0, max_t + 2):
        cum += c.get(t, 0)
        survival.append((total_merged - cum) / total_merged)

    print("\nt : P(tau_early_merge > t | early-merged)")
    for t in [0,1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,70,100,150]:
        if t < len(survival):
            print(f"{t:4d} : {survival[t]:.6f}")

with open('/home/claude/early_merge_times.txt', 'w') as f:
    for t in merge_times:
        f.write(f"{t}\n")

with open('/home/claude/summary.txt', 'w') as f:
    f.write(f"N={N}\nagree={agree}\nearly_merge={early_merge}\n")
