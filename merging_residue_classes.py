def collatz_traj(n, max_steps=400):
    traj = [n]
    while n != 1 and len(traj) < max_steps:
        n = n//2 if n%2==0 else 3*n+1
        traj.append(n)
    return traj

def is_merging_class(r, modulus, num_samples=8, max_check_steps=200):
    merge_steps_seen = set()
    for j in range(2, 2+num_samples):
        n = r + j*modulus
        if n < 2:
            continue
        traj_n = collatz_traj(n, max_check_steps)
        traj_n1 = collatz_traj(n+1, max_check_steps)
        set_n1 = {v:i for i,v in enumerate(traj_n1)}
        merge_step = None
        for idx, v in enumerate(traj_n):
            if v in set_n1:
                merge_step = (idx, set_n1[v])
                break
        if merge_step is None:
            return False, None
        merge_steps_seen.add(merge_step)
    if len(merge_steps_seen) == 1:
        return True, merge_steps_seen.pop()
    return False, merge_steps_seen

print(f"{'modulus':>8} {'merging_residues':>18} {'total_residues':>15} {'merging比率':>12}")
results_by_mod = {}
for k in range(3, 10):
    modulus = 2**k
    merging = []
    for r in range(modulus):   # 全ての余り(偶数も奇数も)をチェック
        ok, info = is_merging_class(r, modulus)
        if ok:
            merging.append((r, info))
    results_by_mod[modulus] = merging
    print(f"{modulus:>8} {len(merging):>18} {modulus:>15} {100*len(merging)/modulus:>11.1f}%")

print("\n=== 各法での「確実に合流するクラス」一覧 (余り, (nから何ステップ, n+1から何ステップ)) ===")
for modulus, merging in results_by_mod.items():
    print(f"mod {modulus} ({len(merging)}個): {merging}")
