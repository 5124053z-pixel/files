import random
from collections import Counter

def collatz_steps(n):
    steps = 0
    while n != 1:
        n = n//2 if n%2==0 else 3*n+1
        steps += 1
    return steps

random.seed(42)

def test(x, L, m_bits, trials=200):
    """xを固定し、ランダムなmに対して steps(2^L*m+x)-steps(m) を計算"""
    diffs = []
    for _ in range(trials):
        m = random.getrandbits(m_bits) | (1 << (m_bits-1))  # ensure exact bit length
        s_m = collatz_steps(m)
        s_big = collatz_steps((m << L) + x)
        diffs.append(s_big - s_m)
    return diffs

test_blocks = [
    ("alternating L=12", 2730, 12),
    ("all-ones L=12", 4095, 12),
    ("random L=12", 2905, 12),
    ("alternating L=20", 699050, 20),
    ("random L=20", 586045, 20),
]

for name, x, L in test_blocks:
    for m_bits in [100, 500]:
        diffs = test(x, L, m_bits, trials=200)
        c = Counter(diffs)
        top_val, top_count = c.most_common(1)[0]
        print(f"{name} (x={x}), ランダムm({m_bits}bit): "
              f"最頻値diff={top_val} (期待L={L}, 一致={top_val==L}), "
              f"頻度={100*top_count/len(diffs):.1f}%, "
              f"平均diff={sum(diffs)/len(diffs):.1f}")
