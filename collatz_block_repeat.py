import random
from collections import Counter
import time

def collatz_steps(n):
    steps = 0
    while n != 1:
        if n % 2 == 0:
            n //= 2
        else:
            n = 3*n + 1
        steps += 1
    return steps

def run_experiment(x, num_iters=80):
    N = x
    prev = collatz_steps(N)
    diffs = []
    for k in range(2, num_iters):
        N = (x << N.bit_length()) | N
        s = collatz_steps(N)
        diffs.append(s - prev)
        prev = s
    return diffs

random.seed(2026)

test_cases = []
# special structured blocks for each length
for L in range(4, 25):
    test_cases.append((f"L={L} all-ones", (1<<L)-1, L))
    test_cases.append((f"L={L} alternating", int('10'*(L//2)+('1' if L%2 else ''), 2), L))
    test_cases.append((f"L={L} single-bit-pair (top+bottom)", (1<<(L-1))|1, L))
    # 3 random odd numbers with top bit set (to guarantee exactly L bits) and bottom bit set (odd)
    for _ in range(2):
        val = (1 << (L-1)) | random.getrandbits(L-2)*2 | 1 if L>=3 else (1<<(L-1))|1
        val |= (1<<(L-1))  # ensure top bit
        val |= 1            # ensure odd
        test_cases.append((f"L={L} random", val, L))

print(f"合計 {len(test_cases)} 個のブロックをテストします\n")

t0 = time.time()
exceptions = []
results_summary = []

for name, x, L in test_cases:
    if x.bit_length() != L:
        continue  # skip malformed
    diffs = run_experiment(x, num_iters=80)
    late = diffs[-30:]  # 定常状態とみなせる後半30個
    c = Counter(late)
    top_val, top_count = c.most_common(1)[0]
    freq = top_count/len(late)
    matches_L = (top_val == L)
    results_summary.append((name, x, L, top_val, freq, matches_L))
    if not matches_L or freq < 0.3:
        exceptions.append((name, x, L, top_val, freq))

elapsed = time.time()-t0
print(f"完了。所要時間 {elapsed:.1f}秒\n")

print(f"{'block':>35} {'L':>4} {'収束先':>8} {'頻度':>7} {'L一致':>6}")
for name, x, L, top_val, freq, matches_L in results_summary:
    print(f"{name:>35} {L:>4} {top_val:>8} {100*freq:>6.1f}% {str(matches_L):>6}")

print(f"\n=== 例外（L以外に収束、または頻度が低いもの） ===")
if exceptions:
    for name, x, L, top_val, freq in exceptions:
        print(f"{name}: x={x}, L={L}, 収束先={top_val} (期待値{L}), 頻度={100*freq:.1f}%")
else:
    print("例外なし。全ケースでL(ブロック自身のビット長)に収束することを確認。")
