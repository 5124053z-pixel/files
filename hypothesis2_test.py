import random
import time
import math

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

def v2(x):
    """2-adic valuation: number of trailing zero bits."""
    if x == 0:
        return 0
    return (x & -x).bit_length() - 1

def trailing_run_length(x):
    """Length of the run of identical bits at the LSB end (0s if x even,
    1s if x odd)."""
    if x & 1 == 0:
        return v2(x)
    # count trailing 1s: invert and count trailing zeros of (~x)&mask
    n = 0
    y = x
    while y & 1 == 1:
        n += 1
        y >>= 1
    return n

BITS = 2048
random.seed(2026)

def random_n(bits):
    return random.getrandbits(bits - 1) | (1 << (bits - 1))

NUM_SAMPLES = 4000
rows = []
t0 = time.time()
attempts = 0
while len(rows) < NUM_SAMPLES:
    attempts += 1
    n = random_n(BITS)
    sX = total_stopping_time(n)
    sY = total_stopping_time(n + 1)
    if sX is None or sY is None or sX != sY:
        continue
    mt = early_merge_time(n, sX)
    if mt is None:
        continue
    rows.append({
        'n_steps': sX,
        'merge_time': mt,
        'v2_n': v2(n),
        'v2_n1': v2(n + 1),
        'run_n': trailing_run_length(n),
        'run_n1': trailing_run_length(n + 1),
        'n_steps_normalized': sX / BITS,
    })

print(f"collected {len(rows)} samples ({attempts} attempts), BITS={BITS}, time={time.time()-t0:.1f}s")

def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs)/n; my = sum(ys)/n
    cov = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    vx = sum((x-mx)**2 for x in xs)
    vy = sum((y-my)**2 for y in ys)
    if vx == 0 or vy == 0:
        return float('nan')
    return cov / math.sqrt(vx*vy)

merge_times = [r['merge_time'] for r in rows]
log_merge_times = [math.log(mt) for mt in merge_times]

for feat in ['n_steps', 'n_steps_normalized', 'v2_n', 'v2_n1', 'run_n', 'run_n1']:
    vals = [r[feat] for r in rows]
    corr_raw = pearson(vals, merge_times)
    corr_log = pearson(vals, log_merge_times)
    print(f"corr({feat:20s}, merge_time)      = {corr_raw:+.4f}")
    print(f"corr({feat:20s}, log(merge_time)) = {corr_log:+.4f}")

# also: mean n_steps_normalized for fast vs slow merge buckets
rows.sort(key=lambda r: r['merge_time'])
n = len(rows)
fast = rows[:n//10]           # bottom 10% (fastest mergers)
slow = rows[-n//10:]          # top 10% (slowest mergers)
mid = rows[n//2 - n//20 : n//2 + n//20]  # median band

def summarize(name, group):
    ns = [r['n_steps_normalized'] for r in group]
    v2ns = [r['v2_n'] for r in group]
    v2n1s = [r['v2_n1'] for r in group]
    print(f"{name}: n={len(group)}  mean(steps/bits)={sum(ns)/len(ns):.4f}  "
          f"mean(v2_n)={sum(v2ns)/len(v2ns):.3f}  mean(v2_n1)={sum(v2n1s)/len(v2n1s):.3f}  "
          f"merge_time range=[{group[0]['merge_time']},{group[-1]['merge_time']}]")

print()
summarize("fastest 10%", fast)
summarize("median band", mid)
summarize("slowest 10%", slow)

with open('/home/claude/hypothesis2_data.csv', 'w') as f:
    f.write("merge_time,n_steps,v2_n,v2_n1,run_n,run_n1\n")
    for r in rows:
        f.write(f"{r['merge_time']},{r['n_steps']},{r['v2_n']},{r['v2_n1']},{r['run_n']},{r['run_n1']}\n")
