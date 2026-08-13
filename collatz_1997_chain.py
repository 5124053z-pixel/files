"""
コラッツ予想 「1997継ぎ足し」探索スクリプト
=====================================================

やっていること:
- 1997からスタートし、「1997」という4桁を延々と継ぎ足していく、たった1本の数列を追いかける
  例(下に継ぎ足す場合): 1997 -> 19971997 -> 199719971997 -> ...
  例(上に継ぎ足す場合): 1997 -> 19971997 -> 199719971997 -> ...(見た目は同じだが桁の伸び方が逆)
- 各回について、コラッツ数列で1に到達するまでのステップ数を記録
- 定期的にcheckpointに保存し、PCを再起動しても続きから再開可能
- バッチ処理 + sleep でCPU負荷を意図的に抑える

止め方: Ctrl+C で安全に停止（次回起動時は続きから再開されます）

注意:
- 桁数がどんどん増えていくので、後半になるほど1回あたりの計算コストが重くなります
  （目安: 前回のベンチマークで1000桁前後は1回あたり数十ミリ秒程度）
- 気になる場合は BATCH_SIZE を減らす、SLEEP_SECONDS を増やす、で調整してください
"""

import os
import sys
import json
import time
import signal

sys.set_int_max_str_digits(0)  # 桁数制限を解除（0=無制限）。巨大な数のstr()変換に必要
_env_runtime = os.environ.get("COLLATZ_MAX_RUNTIME_SECONDS")

try:
    import gmpy2
    _mpz = gmpy2.mpz
    HAS_GMPY2 = True
except ImportError:
    _mpz = int
    HAS_GMPY2 = False

# ==================== 設定 ====================
DIRECTION = "bottom"   # "bottom"=末尾に継ぎ足す(n*10000+1997) / "top"=先頭に継ぎ足す(1997*10^桁+n)
CHECKPOINT_FILE = "collatz_1997_checkpoint.json"
RECORDS_FILE = "collatz_1997_records.log"       # 毎回の結果をすべて記録
ANOMALY_FILE = "collatz_1997_ANOMALY.log"
BATCH_SIZE = 20
SLEEP_SECONDS = 1.0
NICE_LEVEL = 19
# 異常判定の閾値は「桁数 × STEPS_PER_DIGIT_LIMIT」で決める。
# 収束までのステップ数は桁数にほぼ比例して増える（実測: 1桁あたり約23ステップ、
# 最大でも29.7）ため、固定値にすると桁が伸びただけで異常扱いになってしまう。
# 100は実測値の3倍以上の余裕がある。
STEPS_PER_DIGIT_LIMIT = 100
MIN_MAX_STEPS = 10_000   # 桁数が少ないうちの下限
PRINT_EVERY = 100   # コンソールに進捗を表示する頻度（この回数ごとに1回だけ表示。結果は毎回records.logに保存される）
MAX_RUNTIME_SECONDS = None  # Noneなら無制限（PCで動かす用）。GitHub Actionsではワークフロー側から上書きする
if _env_runtime:
    MAX_RUNTIME_SECONDS = float(_env_runtime)
KTABLE_BITS = 18   # kステップテーブルの k。大きいほど速いが構築コストとメモリが増える
# ================================================

running = True


def _build_ktable(k):
    """mod 2^k の値ごとに (mult, add) を求める。
    T^k(n) = (mult*n + add) / 2^k が、n ≡ r (mod 2^k) を満たす全てのnで成り立つ。
    """
    table = {}
    mod = 1 << k
    for r in range(mod):
        n1, n2 = r, r + mod
        for _ in range(k):
            n1 = n1 // 2 if n1 % 2 == 0 else 3 * n1 + 1
            n2 = n2 // 2 if n2 % 2 == 0 else 3 * n2 + 1
        mult = n2 - n1
        add = n1 * mod - mult * r
        table[r] = (mult, add)
    return table


print(f"kステップ高速化テーブルを構築中 (k={KTABLE_BITS}, gmpy2={HAS_GMPY2})...")
_ktable_start = time.time()
_KTABLE = _build_ktable(KTABLE_BITS)
_KMOD = 1 << KTABLE_BITS
print(f"  完了 ({time.time() - _ktable_start:.1f}秒, {len(_KTABLE):,}件)")


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
            p = psutil.Process(os.getpid())
            if sys.platform == "win32":
                p.nice(psutil.IDLE_PRIORITY_CLASS)
                print("プロセス優先度を下げました (IDLE_PRIORITY_CLASS)")
        except ImportError:
            print("優先度調整はスキップしました（psutil未インストール）。動作には影響ありません。")


def append_1997(n):
    """現在のnに『1997』という4桁を継ぎ足した、次の数を返す"""
    if DIRECTION == "bottom":
        return n * 10000 + 1997
    else:  # "top"
        digits = len(str(n))
        return 1997 * (10 ** digits) + n


def step_limit_for(n):
    """nの桁数に応じた異常判定の閾値を返す"""
    return max(MIN_MAX_STEPS, len(str(n)) * STEPS_PER_DIGIT_LIMIT)


def collatz_steps(n, max_steps=None):
    """コラッツ数列が1に到達するまでのステップ数と、途中の最大値を返す。

    n mod 2^k だけで次のkステップ分の変換をまとめて計算する「kステップテーブル」と、
    利用可能ならgmpy2による高速な巨大整数演算を組み合わせている。
    （注: kステップ単位でしか値を見ないため、max_valはジャンプの合間の
    真の最大値を見逃すことがあり、わずかに過小評価される場合がある）
    """
    if max_steps is None:
        max_steps = step_limit_for(n)
    m = _mpz(n)
    steps = 0
    max_val = n
    while m != 1:
        if m >= _KMOD:
            r = int(m & (_KMOD - 1))
            mult, add = _KTABLE[r]
            m = (mult * m + add) >> KTABLE_BITS
            steps += KTABLE_BITS
        else:
            m = m // 2 if m % 2 == 0 else 3 * m + 1
            steps += 1
        m_int = int(m)
        if m_int > max_val:
            max_val = m_int
        if steps >= max_steps:
            return steps, max_val, True
    return steps, max_val, False


def load_checkpoint():
    """チェックポイントを読み込む。

    currentは文字列として保存する。JSONの数値として持つと、読み込み側が
    sys.set_int_max_str_digits()を呼んでいない限り、Python 3.11以降の
    整数⇔文字列変換の桁数制限（既定4300桁）に引っかかって例外になるため。
    （旧形式のintで保存されたファイルもそのまま読める）
    """
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            data = json.load(f)
        data["current"] = int(data["current"])
        print(f"チェックポイントを読み込みました: {data['iteration']}回目, "
              f"n={str(data['current'])[:20]}...({len(str(data['current']))}桁)")
        return data
    return {
        "current": 1997,
        "iteration": 0,
    }


def save_checkpoint(data):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"current": str(data["current"]), "iteration": data["iteration"]}, f)


def log_result(iteration, n, steps, max_val):
    with open(RECORDS_FILE, "a") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | "
            f"iteration={iteration} | digits={len(str(n))} | "
            f"steps={steps} | max_val_digits={len(str(max_val))}\n"
        )


def log_anomaly(n, steps):
    digits = len(str(n))
    with open(ANOMALY_FILE, "a") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | "
            f"digits={digits} | limit={step_limit_for(n):,} | "
            f"n={n} が {steps:,}ステップを超えても1に到達しませんでした。"
            f"ループまたは発散の可能性があります。要確認。\n"
        )


def main():
    lower_priority()
    state = load_checkpoint()
    n = state["current"]
    iteration = state["iteration"]
    start_time = time.time()

    print(f"探索を開始します（継ぎ足し方向: {DIRECTION}）。停止するには Ctrl+C を押してください。")
    if MAX_RUNTIME_SECONDS:
        print(f"制限時間: {MAX_RUNTIME_SECONDS/3600:.1f}時間で自動的に安全終了します。")
    if len(str(n)) > 30:
        print(f"開始位置: {iteration}回目, n={str(n)[:30]}...({len(str(n))}桁)")
    else:
        print(f"開始位置: {iteration}回目, n={n}")

    global running
    while running:
        if MAX_RUNTIME_SECONDS and time.time() - start_time >= MAX_RUNTIME_SECONDS:
            print(f"制限時間({MAX_RUNTIME_SECONDS/3600:.1f}時間)に達したため、安全に終了します。")
            running = False
            break

        for _ in range(BATCH_SIZE):
            if not running:
                break

            steps, max_val, hit_limit = collatz_steps(n)

            if hit_limit:
                log_anomaly(n, steps)
                print("=" * 60)
                print(f"⚠️  異常検知: n={str(n)[:30]}... が {steps:,}ステップを超えても収束しません")
                print(f"    {ANOMALY_FILE} に記録し、安全のため探索を停止します。")
                print("=" * 60)
                state["current"] = n
                state["iteration"] = iteration
                save_checkpoint(state)
                return

            iteration += 1
            log_result(iteration, n, steps, max_val)

            if iteration % PRINT_EVERY == 0:
                digits = len(str(n))
                n_display = str(n) if digits <= 20 else f"{str(n)[:20]}..."
                print(f"{iteration:4d}回目 ({digits:5d}桁, n={n_display}) -> {steps:,}ステップ")

            n = append_1997(n)

        state["current"] = n
        state["iteration"] = iteration
        save_checkpoint(state)

        time.sleep(SLEEP_SECONDS)

    save_checkpoint(state)
    print(f"終了しました。次回は{iteration}回目({len(str(n))}桁)から再開します。")


if __name__ == "__main__":
    main()
