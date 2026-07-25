import random
import time

def early_merge_time(n, cap):
    x, y = n, n + 1
    for t in range(1, cap):
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        y = y >> 1 if y & 1 == 0 else 3 * y + 1
        if x == y:
            return t
    return None

def total_stopping_time(n, hard_cap=2_000_000):
    t = 0
    x = n
    while x != 1:
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        t += 1
        if t > hard_cap:
            return None
    return t

BITS = 4096  # deliberately large, to make clear any small K found is NOT
             # just "K happens to be close to the full bit length"
NUM_CONTROLS_PER_TEST = 4
MAX_K = 60

random.seed(777)

def random_n(bits):
    return random.getrandbits(bits - 1) | (1 << (bits - 1))

def explained_by_modulus(n, K, mt):
    """Does fixing n's low K bits (with fresh random high bits) reliably
    reproduce the same merge time mt, across NUM_CONTROLS_PER_TEST trials?"""
    low_bits = n % (1 << K)
    for _ in range(NUM_CONTROLS_PER_TEST):
        high = random.getrandbits(BITS - K)
        m = (high << K) | low_bits
        m |= (1 << (BITS - 1))
        sm = total_stopping_time(m)
        sm1 = total_stopping_time(m + 1)
        if sm is None or sm1 is None or sm != sm1:
            return False
        mtm = early_merge_time(m, sm)
        if mtm != mt:
            return False
    return True

def minimal_explaining_K(n, mt):
    """Binary search smallest K in [3, MAX_K] such that low K bits explain
    the merge. Returns K, or None if not found within MAX_K."""
    lo, hi = 3, MAX_K
    if not explained_by_modulus(n, hi, mt):
        return None  # not explained even by MAX_K bits -- true "slow" case
    while lo < hi:
        mid = (lo + hi) // 2
        if explained_by_modulus(n, mid, mt):
            hi = mid
        else:
            lo = mid + 1
    return lo

# collect "fast" samples (merge_time < 20) and find minimal explaining K
NUM_SAMPLES = 150
collected = []
attempts = 0
t0 = time.time()
while len(collected) < NUM_SAMPLES:
    attempts += 1
    n = random_n(BITS)
    s = total_stopping_time(n)
    s1 = total_stopping_time(n + 1)
    if s is None or s1 is None or s != s1:
        continue
    mt = early_merge_time(n, s)
    if mt is None or mt >= 20:
        continue
    collected.append((n, mt))

print(f"collected {len(collected)} fast (t<20) samples from {attempts} attempts at BITS={BITS}, time={time.time()-t0:.1f}s")

Ks = []
t1 = time.time()
for n, mt in collected:
    K = minimal_explaining_K(n, mt)
    Ks.append((mt, K))

print(f"minimal-K search done, time={time.time()-t1:.1f}s")

found = [K for mt, K in Ks if K is not None]
notfound = [mt for mt, K in Ks if K is None]

print(f"\nminimal explaining K found for {len(found)}/{len(Ks)} samples")
if found:
    print(f"min K: {min(found)}, max K: {max(found)}, mean K: {sum(found)/len(found):.2f}")
    from collections import Counter
    c = Counter(found)
    print("distribution of minimal K:")
    for k in sorted(c):
        print(f"  K={k:3d}: {c[k]} samples")
if notfound:
    print(f"\n{len(notfound)} samples not explained even by K={MAX_K} (merge_times: {notfound})")

with open('/home/claude/minimal_K.csv', 'w') as f:
    f.write("merge_time,minimal_K\n")
    for mt, K in Ks:
        f.write(f"{mt},{K if K is not None else -1}\n")
