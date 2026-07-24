def collatz_traj(n, max_steps=400):
    traj = [n]
    while n != 1 and len(traj) < max_steps:
        n = n//2 if n%2==0 else 3*n+1
        traj.append(n)
    return traj

def is_merging_class_general(x, L, r, modulus, num_samples=6, max_check_steps=300):
    """m ≡ r (mod modulus) のとき、mと(2^L*m+x)が常に同じステップ数で合流するか"""
    merge_steps_seen = set()
    for j in range(2, 2+num_samples):
        m = r + j*modulus
        if m < 1:
            continue
        N = (m << L) + x
        traj_m = collatz_traj(m, max_check_steps)
        traj_N = collatz_traj(N, max_check_steps)
        set_N = {v:i for i,v in enumerate(traj_N)}
        merge = None
        for idx, v in enumerate(traj_m):
            if v in set_N:
                merge = (idx, set_N[v])
                break
        if merge is None:
            return False
        merge_steps_seen.add(merge)
    return len(merge_steps_seen) == 1

def merging_fraction(x, L, k_range):
    results = []
    for k in k_range:
        modulus = 2**k
        count = sum(1 for r in range(modulus) if is_merging_class_general(x, L, r, modulus))
        frac = count/modulus
        results.append((modulus, count, frac))
    return results

print("=== alternating x (L=12, x=2730) の合流クラス率 ===")
for modulus, count, frac in merging_fraction(2730, 12, range(1,7)):
    print(f"  mod {modulus}: {count}/{modulus} = {100*frac:.1f}%")

print("\n=== ランダムっぽい x (L=12, x=2905) の合流クラス率 ===")
for modulus, count, frac in merging_fraction(2905, 12, range(1,7)):
    print(f"  mod {modulus}: {count}/{modulus} = {100*frac:.1f}%")

print("\n=== ランダムx をさらに大きいmodulusまで調べる ===")
for modulus, count, frac in merging_fraction(2905, 12, range(7,10)):
    print(f"  mod {modulus}: {count}/{modulus} = {100*frac:.2f}%")
