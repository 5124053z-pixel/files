from collections import Counter

def collatz_steps(n):
    steps = 0
    while n != 1:
        if n % 2 == 0:
            n //= 2
        else:
            n = 3*n + 1
        steps += 1
    return steps

def run_experiment(x, num_iters=70):
    N = x
    prev = collatz_steps(N)
    diffs = []
    for k in range(2, num_iters):
        N = (x << N.bit_length()) | N
        s = collatz_steps(N)
        diffs.append(s - prev)
        prev = s
    return diffs

def late_freq_of_L(x, L, num_iters=70, tail=30):
    diffs = run_experiment(x, num_iters)
    late = diffs[-tail:]
    c = Counter(late)
    return c.get(L, 0) / len(late), c.most_common(3)

# 「数ステップ以内に、ぴったり2のべき乗に着地する数」を逆算して探す
# power_of_2の逆コラッツをたどって、いくつかステップ手前の値を求める
def find_preimages_landing_on_power_of_2(target_power, steps_back, max_results=5):
    """2^target_power から steps_back ステップ逆算して、その手前の数を集める"""
    current_set = {1 << target_power}
    for _ in range(steps_back):
        new_set = set()
        for v in current_set:
            # 偶数への逆操作: v -> 2v (常に可能)
            new_set.add(2*v)
            # 奇数への逆操作: v = 3u+1 -> u = (v-1)/3 (v≡1 mod3 かつ (v-1)/3が奇数の場合のみ)
            if v % 3 == 1:
                u = (v-1)//3
                if u % 2 == 1 and u > 0:
                    new_set.add(u)
        current_set = new_set
    return list(current_set)[:max_results]

results = []
print("=== 「数ステップでぴったり2のべき乗に着地する」数を探索して検証 ===\n")
for target_power in [12, 14, 16, 18]:
    for steps_back in [2, 3, 4]:
        candidates = find_preimages_landing_on_power_of_2(target_power, steps_back, max_results=2)
        for x in candidates:
            if x <= 1:
                continue
            L = x.bit_length()
            if L < 6:
                continue
            own_steps = collatz_steps(x)
            freq, top3 = late_freq_of_L(x, L)
            results.append((x, L, own_steps, steps_back, target_power, freq))
            print(f"x={x:>10} (L={L:>2}), 目標2^{target_power}まで{steps_back}ステップ, own_steps={own_steps:>4}, "
                  f"freq(diff==L)={100*freq:>5.1f}%  top3={top3}")

print(f"\n全{len(results)}件の平均freq: {100*sum(r[5] for r in results)/len(results):.1f}%")
