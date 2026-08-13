"""
コラッツ予想 ループ実在探索スクリプト（存在側を攻める）
=====================================================

collatz_cycle_search.py との役割分担:
- あちらは「ループが存在するなら (k,S) はこの18組のどれか」まで絞る。
  ただし非存在の証明は Baker の理論が要るため、そこで手が止まる。
- こちらは非存在の証明を諦め、「あるなら見つける」側だけを狙う。
  非存在の証明は無限個を相手にするので有限回では終わらないが、
  存在の発見は1個見つければ即座に決着がつく。この非対称性を使う。

素朴にやると破綻する:
  指数列 (e_0..e_{k-1}) を総当たりすると k=306 で 10^137 通り。
  さらに n_0 = C/D が整数になる確率が ~1/D ~ 10^-143 しかない。
  「指数列を選んで n_0 を計算する」方向は完全に無理。

そこで向きを逆にする:
  n_0 を先に決めて実際に加速コラッツ写像を回せば、割り算は必ず割り切れる
  （e_i は「2で割れるだけ割る」と定義されるため）。整除条件が自動的に満たされる。
  あとは k ステップ後に出発点に戻るかを見るだけ。戻れば本物のループ。

  つまり「n_0 を1つ選ぶ → k回回す → 戻ったか？」の繰り返し。
  1回の判定は k ステップの計算だけで済み、整除の奇跡を待つ必要がない。

それでも残る問題（正直に）:
  候補の (k,S) から逆算すると、ループが存在する場合の n_0 は
      k=306   -> 約149桁
      k=665   -> 約322桁
      k=15601 -> 約7448桁
  という大きさになる。149桁の数は 10^149 個あり、そこから正解を引く確率は
  ランダム探索と同じく事実上ゼロ。
  つまり本スクリプトも「当たれば大発見、当たらない可能性が圧倒的」という
  宝くじであることに変わりはない。ランダム探索との違いは、
    - 対象が「1に収束するか」ではなく「ループを成すか」という強い条件
    - k と S が候補に絞られているので、回すステップ数が確定している
    - 走査を構造的に行える（下記）
  という点だけ。

枝刈り:
  n_0 を「ループ内の最小要素」として探すので、軌道が一度でも n_0 を
  下回ったらループになりえず、その場で打ち切れる。加速写像は1ステップ
  あたり平均 log10(3/4) ≈ -0.125 桁ずつ縮むので、ほとんどの n_0 は
  数ステップで脱落し、k回まわしきる必要がない。

  なお当初は「2^71 以下に落ちたら失格」という枝刈りを入れていたが、
  これは一度も発動しなかった。149桁を22桁まで落とすには約1000ステップ
  かかるのに、k=306 では306ステップしか回さないため。
  n_0 との比較に変えたことで、枝刈りとして実際に機能するようになった。

止め方: Ctrl+C で安全に停止（次回起動時は続きから再開）
"""

import os
import sys
import json
import time
import math
import random
import signal

sys.set_int_max_str_digits(0)

try:
    import gmpy2
    _mpz = gmpy2.mpz
    HAS_GMPY2 = True
except ImportError:
    _mpz = int
    HAS_GMPY2 = False

# ==================== 設定 ====================
CHECKPOINT_FILE = "collatz_hunt_checkpoint.json"
RECORDS_FILE = "collatz_hunt_records.log"
FOUND_FILE = "collatz_hunt_FOUND.log"   # 本物のループが見つかったらここに書く

VERIFIED_BOUND = 2 ** 71    # 総当たり検証済み。これ以下に落ちたら失格
BATCH_SIZE = 200            # 1バッチで試す n_0 の個数
SLEEP_SECONDS = 1.0
NICE_LEVEL = 19
PRINT_EVERY = 20            # バッチ何回ごとに進捗を出すか
_env_runtime = os.environ.get("COLLATZ_MAX_RUNTIME_SECONDS")
MAX_RUNTIME_SECONDS = float(_env_runtime) if _env_runtime else None

# 探索対象の (k, S)。collatz_cycle_search.py が出した候補のうち
# 現実的に回せる小さいものから。k が大きいほど1回の判定が重くなる。
# S は ceil(k*log2(3)) 以上でなければ 2^S > 3^k を満たさず、n_0 が正にならない。
TARGETS = [
    (306, 485),
    (665, 1055),
]
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


def expected_n0_digits(k, S):
    """(k,S) でループがある場合の n_0 の桁数の目安を返す。

    n_0 = C/D、C ≈ (6/5)·3^{k-1}·2^{S-1}、D = 2^S - 3^k から計算する。
    """
    D = (1 << S) - 3 ** k
    if D <= 0:
        return None
    n0 = (6 * 3 ** (k - 1) * (1 << (S - 1))) // (5 * D)
    return len(str(n0))


def try_cycle(n0, k):
    """n_0 から加速コラッツ写像を k 回回し、出発点に戻るか調べる。

    戻り値 (is_cycle, S_actual, steps_done)
      is_cycle   : k回でちょうど n_0 に戻ったか
      S_actual   : その間に 2 で割った回数の合計
      steps_done : 失格するまでに進めたステップ数（k なら最後まで到達）

    ここでは指数列を作らず実際に写像を適用するので、割り算は必ず割り切れる。
    整除条件を人為的に満たす必要がないのが、この向きの利点。

    枝刈り: n_0 を「ループ内の最小要素」として探すので、途中で n_0 を
    下回ったらその時点で失格にできる。2^71 との比較では k=306 のとき
    149桁を22桁まで落とすのに約1000ステップ必要で306回では届かず、
    一度も発動しなかった。n_0 との比較なら平均数ステップで脱落するため、
    枝刈りとして実際に機能する。
    """
    n = _mpz(n0)
    target = _mpz(n0)
    S_actual = 0
    for i in range(k):
        n = 3 * n + 1
        while n % 2 == 0:
            n //= 2
            S_actual += 1
        if n == target:
            # k回より早く戻った = 周期がkの約数。k回でも戻るので成立
            return True, S_actual, i + 1
        # ループの最小要素を n_0 と仮定しているので、下回ったら失格
        if n < target:
            return False, S_actual, i + 1
    return False, S_actual, k


def random_n0(digits):
    """指定桁数のランダムな奇数を返す"""
    lo = 10 ** (digits - 1)
    hi = 10 ** digits - 1
    return random.randint(lo, hi) | 1


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        print(f"チェックポイントを読み込みました: "
              f"{data['tested']:,}個テスト済み, 最長{data['best_return']}ステップ生存")
        return data
    return {"tested": 0, "best_return": 0, "target_idx": 0,
            "fell_below": 0, "survived": 0}


def save_checkpoint(data):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f)


def log_found(n0, k, S_actual):
    """本物のループを発見した場合（これが出たらコラッツ予想の反例）"""
    with open(FOUND_FILE, "a", encoding="utf-8") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | "
            f"★★★ ループ発見 ★★★ k={k} | S={S_actual} | "
            f"n_0={n0}\n"
        )


def log_progress(state, k, S, digits, elapsed):
    with open(RECORDS_FILE, "a", encoding="utf-8") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | k={k} | S={S} | "
            f"n0桁数={digits} | テスト数={state['tested']} | "
            f"脱落={state['fell_below']} | 完走={state['survived']} | "
            f"{elapsed:.0f}秒\n"
        )


def main():
    lower_priority()
    state = load_checkpoint()
    start_time = time.time()

    print("=" * 72)
    print("ループ実在探索（存在側のみを狙う）")
    print(f"  既検証境界: 2^71 = {VERIFIED_BOUND:,}（これ以下に落ちたら失格）")
    print(f"  gmpy2: {HAS_GMPY2}")
    print("=" * 72)
    for k, S in TARGETS:
        d = expected_n0_digits(k, S)
        print(f"  対象 k={k}, S={S} → 探すべき n_0 は約 {d} 桁 "
              f"（その空間は約 10^{d} 個）")
    print("=" * 72)
    print("注意: 空間が巨大なため、当たる見込みは事実上ありません。")
    print("      『あるなら見つかる』方式であり、見つからないことは")
    print("      『ない』ことの証明にはなりません。")
    print("=" * 72)
    if MAX_RUNTIME_SECONDS:
        print(f"制限時間: {MAX_RUNTIME_SECONDS/3600:.2f}時間")
    print()

    global running
    batch_count = 0
    k, S = TARGETS[state["target_idx"] % len(TARGETS)]
    digits = expected_n0_digits(k, S)

    while running:
        if MAX_RUNTIME_SECONDS and time.time() - start_time >= MAX_RUNTIME_SECONDS:
            print("制限時間に達したため、安全に終了します。")
            break

        for _ in range(BATCH_SIZE):
            if not running:
                break
            n0 = random_n0(digits)
            is_cycle, S_actual, steps_done = try_cycle(n0, k)
            state["tested"] += 1
            if steps_done > state["best_return"]:
                state["best_return"] = steps_done

            if is_cycle:
                # ここに来たらコラッツ予想の反例。最重要。
                log_found(n0, k, S_actual)
                print("=" * 72)
                print("★★★ ループを発見しました ★★★")
                print(f"  k={k}, S={S_actual}")
                print(f"  n_0={n0}")
                print(f"  {FOUND_FILE} に記録しました。")
                print("=" * 72)
                save_checkpoint(state)
                return

            if steps_done < k:
                # 途中で n_0 を下回った = ループになりえない
                state["fell_below"] += 1
            else:
                # k回まわっても一度も n_0 を下回らなかった = 珍しい
                state["survived"] += 1

        batch_count += 1
        save_checkpoint(state)

        if batch_count % PRINT_EVERY == 0:
            elapsed = time.time() - start_time
            rate = state["tested"] / elapsed if elapsed else 0
            print(f"k={k} | テスト {state['tested']:,}個 | "
                  f"n0未満で脱落 {state['fell_below']:,} | "
                  f"k回完走 {state['survived']:,} | "
                  f"{rate:.0f}個/秒 | {elapsed:.0f}秒")
            log_progress(state, k, S, digits, elapsed)

        time.sleep(SLEEP_SECONDS)

    save_checkpoint(state)
    elapsed = time.time() - start_time
    print()
    print("=" * 72)
    print(f"終了しました。")
    print(f"  テストした n_0: {state['tested']:,} 個")
    print(f"  途中で n_0 を下回り脱落: {state['fell_below']:,} 個")
    print(f"  k回完走した: {state['survived']:,} 個")
    print(f"  ループ発見: 0 件")
    print(f"  （探索空間 10^{digits} に対し {state['tested']:,} 個 = "
          f"割合 10^-{digits - len(str(state['tested'])) if state['tested'] else digits} 程度）")
    print("=" * 72)


if __name__ == "__main__":
    main()
