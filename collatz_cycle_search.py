"""
コラッツ予想 ループ（巡回）構造探索スクリプト
=====================================================

これまでのランダム探索との違い:
- ランダム探索は「適当な巨大数を選んで1に落ちるか見る」方式だった。
  反例の下限を押し上げず、10^1026 の空間から数千万個引いても割合は 10^-1019。
- こちらは「ループが存在するとしたら (k,S) はどうでなければならないか」
  という制約から候補を直接構成する。候補は無限ではなく、数個しかない。

理論:
  奇数 n に対する加速コラッツ写像を n -> (3n+1)/2^e （e は割り切れる限り最大）
  とする。奇数を k 個含むループがあり、i番目の割り算の指数を e_i、S = Σe_i と
  すると、

      n_0 = C / (2^S - 3^k),   C = Σ_{i=0}^{k-1} 3^{k-1-i} * 2^{S - E_i}
      （E_i = e_0 + ... + e_i）

  n_0 は正の奇数でなければならないので 2^S > 3^k、すなわち S/k > log2(3)。

絞り込みの核心:
  n_0 = C/D の D = 2^S - 3^k が大きいほど n_0 は小さくなり、総当たり検証済み
  領域（2^71 まで）に落ちて反例になりえない。つまり D が「異常に小さい」
  (k,S) だけが候補として残る。

  D が小さい ⇔ S/k が log2(3) に極めて近い ⇔ S/k が log2(3) の連分数
  収束分数（またはその近く）である。収束分数は疎にしか現れないので、
  候補は k <= 5,000,000 でも十数個しかない。

  ただし正直に言うと、本スクリプトが使う n_0 の上界
      n_0 <= (6/5) * 3^{k-1} * 2^{S-1} / D
  は緩すぎて、実際にはどの候補も棄却できない（n_0上限が149桁〜9万桁になり、
  2^71 = 22桁を大きく超えるため）。候補を数個に絞るところまでは正しく効くが、
  そこから先を潰すには Baker の対数一次形式の理論による下界評価が必要で、
  それが Steiner / Simons-de Weger / Hercher の仕事にあたる。
  本スクリプトは「候補生成までは自力で到達でき、その先には別の道具が要る」
  という境界を実際に示すもの。

  候補となる k は k <= 5,000,000 の範囲で
      k = 2, 5, 12, 41, 53, 306, 665, 15601, 31867, 79335, 111202, 190537, ...
  のわずか十数個。k <= 91 が排除済みなので、実際に残るのは k = 306 以降。
  ランダム探索が 10^1026 の空間を相手にしていたのに対し、こちらは十数個。
  これが「構造を使う」ということ。

正直な注記:
  この方向は既に専門家が徹底的にやっている（Steiner 1977, Simons 1990,
  Simons-de Weger 2005, Hercher 2023）。本スクリプトが新しい定理を出す見込みは
  ない。目的は「構造を使うと空間がどれだけ縮むか」を実際に確かめること。
  ランダム探索よりは筋が良い、という位置づけ。

止め方: Ctrl+C で安全に停止（次回起動時は続きから再開）
"""

import os
import sys
import json
import time
import math
import signal

sys.set_int_max_str_digits(0)

try:
    import gmpy2
    HAS_GMPY2 = True
except ImportError:
    HAS_GMPY2 = False

# ==================== 設定 ====================
CHECKPOINT_FILE = "collatz_cycle_checkpoint.json"
RECORDS_FILE = "collatz_cycle_records.log"
ANOMALY_FILE = "collatz_cycle_ANOMALY.log"   # 本物のループ候補が出たらここに書く

K_MIN = 92              # Hercher 2023 が k<=91 を排除済み。その次から
K_MAX = 5_000_000       # 候補を探すループ長の上限
VERIFIED_BOUND = 2 ** 71    # 総当たり検証済みの下限
SEARCH_RADIUS = 2       # 各収束分数 S/k の周りで S を ±この範囲だけ試す
SLEEP_SECONDS = 1.0
NICE_LEVEL = 19
_env_runtime = os.environ.get("COLLATZ_MAX_RUNTIME_SECONDS")
MAX_RUNTIME_SECONDS = float(_env_runtime) if _env_runtime else None
# ================================================

running = True
LOG2_3 = math.log2(3)


def handle_stop(signum, frame):
    global running
    print("\n停止シグナルを受け取りました。安全に終了します...")
    running = False


signal.signal(signal.SIGINT, handle_stop)
signal.signal(signal.SIGTERM, handle_stop)


def lower_priority():
    try:
        os.nice(NICE_LEVEL)
        print(f"プロセス優先度を下げました (nice={NICE_LEVEL})")
    except (AttributeError, OSError):
        try:
            import psutil
            if sys.platform == "win32":
                psutil.Process(os.getpid()).nice(psutil.IDLE_PRIORITY_CLASS)
                print("プロセス優先度を下げました (IDLE_PRIORITY_CLASS)")
        except ImportError:
            print("優先度調整はスキップしました（psutil未インストール）。")


def convergents(limit_k):
    """log2(3) の連分数収束分数 (k, S) を、k <= limit_k の範囲で返す。

    S/k が log2(3) に近いほど D = 2^S - 3^k が小さくなり、ループ候補になる。
    収束分数はその「最も近い有理数」なので、ここだけ調べればよい。
    """
    x = LOG2_3
    terms = []
    v = x
    for _ in range(60):
        i = int(v)
        terms.append(i)
        frac = v - i
        if frac < 1e-15:
            break
        v = 1.0 / frac

    out = []
    h1, h0 = 1, 0   # 分子 S
    k1, k0 = 0, 1   # 分母 k
    for a in terms:
        h1, h0 = a * h1 + h0, h1
        k1, k0 = a * k1 + k0, k1
        if k1 > limit_k:
            break
        if k1 >= 1:
            out.append((k1, h1))
    return out


def analyze_candidate(k, S):
    """(k,S) について、ループが存在しうるかを評価する。

    戻り値 dict:
      D        : 2^S - 3^k （小さいほど危険＝候補として生き残る）
      n0_upper : ループ内最小要素 n_0 の上限（緩い評価）
      rejected : 既検証境界により棄却できたか
      reason   : 棄却理由
    """
    if S <= 0 or k <= 0:
        return None
    # 2^S > 3^k でなければ n_0 が正にならない
    pow2 = 1 << S
    pow3 = 3 ** k
    D = pow2 - pow3
    if D <= 0:
        return {"k": k, "S": S, "D": D, "n0_upper": None,
                "rejected": True, "reason": "2^S <= 3^k （n_0が正にならない）"}

    # n_0 の上限を求める。
    #   C = Σ_{i} 3^{k-1-i} 2^{S-E_i} は E_i を最小（e_i を後ろに寄せる）に
    #   したとき最大になり、そのとき term_i = 3^{k-1-i} 2^{S-i-1}。
    #   比 term_{i+1}/term_i = 1/6 なので等比級数となり、和は支配項 term_0 の
    #   6/5 倍に収束する（数値実験でも全kで C/term_0 = 1.2 を確認）。
    #     C <= (6/5) * 3^{k-1} * 2^{S-1}
    #   よって n_0 = C/D <= (6/5) * 3^{k-1} * 2^{S-1} / D。
    n0_upper = (6 * 3 ** (k - 1) * (1 << (S - 1))) // (5 * D)

    if n0_upper <= VERIFIED_BOUND:
        return {"k": k, "S": S, "D": D, "n0_upper": n0_upper,
                "rejected": True,
                "reason": f"n_0上限 {n0_upper} <= 2^71（総当たり検証済み領域）"}

    return {"k": k, "S": S, "D": D, "n0_upper": n0_upper,
            "rejected": False, "reason": "この上界では棄却できず"}


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        print(f"チェックポイントを読み込みました: 候補{data['idx']}個目から再開")
        return data
    return {"idx": 0, "checked": 0, "survivors": 0}


def save_checkpoint(data):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f)


def log_result(r, is_survivor):
    path = ANOMALY_FILE if is_survivor else RECORDS_FILE
    with open(path, "a", encoding="utf-8") as f:
        n0d = len(str(r["n0_upper"])) if r["n0_upper"] is not None else 0
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | k={r['k']} | S={r['S']} | "
            f"gap={r['S']/r['k'] - LOG2_3:.6e} | "
            f"D桁数={len(str(abs(r['D'])))} | n0上限桁数={n0d} | "
            f"{'★未棄却' if is_survivor else '棄却'}: {r['reason']}\n"
        )


def main():
    lower_priority()
    state = load_checkpoint()
    start_time = time.time()

    print("=" * 72)
    print("ループ構造探索（連分数による候補生成）")
    print(f"  k <= 91 は Hercher 2023 により排除済み → k >= {K_MIN} を対象")
    print(f"  既検証境界: 2^71 = {VERIFIED_BOUND:,}")
    print(f"  gmpy2: {HAS_GMPY2}")
    print("=" * 72)

    # 候補 (k,S) を生成: 収束分数とその近傍
    conv = convergents(K_MAX)
    candidates = []
    seen = set()
    for k, S in conv:
        if k < K_MIN:
            continue
        for ds in range(-SEARCH_RADIUS, SEARCH_RADIUS + 1):
            S2 = S + ds
            # S/k <= log2(3) では 2^S <= 3^k となり n_0 が正にならない。
            # ループの候補ですらないので、ここで除外しておく
            # （残さないと「棄却」件数がこの自明ケースで埋まってしまう）
            if S2 <= k * LOG2_3:
                continue
            key = (k, S2)
            if key in seen:
                continue
            seen.add(key)
            candidates.append(key)
    candidates.sort()

    print(f"生成された候補 (k,S): {len(candidates)} 件")
    print(f"  （ランダム探索が相手にしていた空間: 約 10^1026 個）")
    print()
    print(f"{'k':>10} {'S':>10} {'S/k-log2(3)':>14} {'D桁数':>8} {'判定':>10}")
    print("-" * 72)

    global running
    for i in range(state["idx"], len(candidates)):
        if not running:
            break
        if MAX_RUNTIME_SECONDS and time.time() - start_time >= MAX_RUNTIME_SECONDS:
            print("\n制限時間に達したため、安全に終了します。")
            break

        k, S = candidates[i]
        r = analyze_candidate(k, S)
        state["checked"] += 1
        if r is None:
            state["idx"] = i + 1
            continue

        is_survivor = not r["rejected"]
        if is_survivor:
            state["survivors"] += 1
        log_result(r, is_survivor)

        gap = S / k - LOG2_3
        mark = "★未棄却" if is_survivor else "棄却"
        print(f"{k:>10} {S:>10} {gap:>14.3e} {len(str(abs(r['D']))):>8} {mark:>10}")

        state["idx"] = i + 1
        save_checkpoint(state)

    save_checkpoint(state)
    print("-" * 72)
    print(f"検査した候補: {state['checked']} 件")
    print(f"棄却できなかった候補: {state['survivors']} 件")
    if state["survivors"]:
        print(f"  → {ANOMALY_FILE} を確認してください")
        print("     （「ループが見つかった」ではなく「素朴な議論では潰せなかった」の意味）")
    print("=" * 72)


if __name__ == "__main__":
    main()
