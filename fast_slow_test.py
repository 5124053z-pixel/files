import random
import time

def early_merge_time(n, cap):
    """First t in [1, cap) with T^t(n) == T^t(n+1), raw map."""
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

BITS = 1024
K = 24  # modulus 2^K used to test "small-modulus explains it" hypothesis
MOD = 1 << K
NUM_SAMPLES = 300
NUM_CONTROLS = 3  # controls per sample, sharing low K bits, different high bits

random.seed(12345)

def random_n(bits):
    return random.getrandbits(bits - 1) | (1 << (bits - 1))

results = []  # (merge_time, explained_bool, num_controls_matching)
t0 = time.time()

collected = 0
attempts = 0
while collected < NUM_SAMPLES:
    attempts += 1
    n = random_n(BITS)
    s = total_stopping_time(n)
    s1 = total_stopping_time(n + 1)
    if s is None or s1 is None or s != s1:
        continue
    mt = early_merge_time(n, s)
    if mt is None:
        continue  # shouldn't happen per prior findings, but guard anyway

    low_bits = n % MOD

    # build controls: same low K bits, fresh random high bits
    matches = 0
    for _ in range(NUM_CONTROLS):
        high = random.getrandbits(BITS - K)
        m = (high << K) | low_bits
        # ensure full bit length roughly preserved (does not touch low K bits)
        m |= (1 << (BITS - 1))
        sm = total_stopping_time(m)
        sm1 = total_stopping_time(m + 1)
        if sm is not None and sm1 is not None and sm == sm1:
            mtm = early_merge_time(m, sm)
            if mtm == mt:
                matches += 1

    explained = (matches == NUM_CONTROLS)
    results.append((mt, explained, matches))
    collected += 1

print(f"collected {collected} samples in {attempts} attempts, time={time.time()-t0:.1f}s")

# bucket by merge_time and report "explained by mod 2^K" fraction per bucket
results.sort(key=lambda r: r[0])
buckets = [(0,10),(10,20),(20,40),(40,80),(80,160),(160,400),(400,100000)]
print(f"\nmodulus tested: 2^{K}")
print(f"{'t range':>15} {'n':>5} {'frac fully explained (all controls match)':>42} {'avg matches/{}'.format(NUM_CONTROLS):>18}")
for lo, hi in buckets:
    sub = [r for r in results if lo <= r[0] < hi]
    if not sub:
        continue
    frac_explained = sum(1 for r in sub if r[1]) / len(sub)
    avg_matches = sum(r[2] for r in sub) / len(sub)
    print(f"{lo:>6}-{hi:<7} {len(sub):>5} {frac_explained:>42.3f} {avg_matches:>18.3f}")

with open('/home/claude/fast_slow_explain.csv', 'w') as f:
    f.write("merge_time,explained,matches\n")
    for mt, exp, m in results:
        f.write(f"{mt},{int(exp)},{m}\n")
