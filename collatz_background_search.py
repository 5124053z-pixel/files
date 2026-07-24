"""
コラッツ予想 バックグラウンド探索スクリプト（低負荷版）
=====================================================

やっていること:
- 毎回ランダムな大きな奇数を選び、コラッツ数列を計算して1に収束することを確認
- ステップ数（軌道の長さ）が過去最高を更新したら records.log に記録
- 定期的に進捗を checkpoint.json に保存し、PCを再起動しても続きから再開可能
- バッチ処理 + sleep でCPU負荷を意図的に抑える（他の作業を邪魔しない）

止め方: Ctrl+C で安全に停止（次回起動時は続きから再開されます）

負荷調整:
- BATCH_SIZE を小さくする / SLEEP_SECONDS を大きくすると、さらに軽くなります
- os.nice() で優先度も下げているので、他のアプリが優先されます
"""

import os
import sys
import json
import time
import signal
import random

sys.set_int_max_str_digits(0)  # 桁数制限を解除（0=無制限）。巨大な数のstr()変換に必要

try:
    import gmpy2
    _mpz = gmpy2.mpz
    HAS_GMPY2 = True
except ImportError:
    _mpz = int
    HAS_GMPY2 = False

# ==================== 設定（ここを調整して負荷をコントロール） ====================
START_NUMBER = 1997 * 10 ** 1022  # 1.997 * 10^1025（初回のみ使用。以降はcheckpointから再開）
BATCH_SIZE = 30                 # 1回の作業でチェックする個数（桁数が大きいので少なめに設定）
SLEEP_SECONDS = 1.0             # バッチ間の休憩時間（長いほど軽い）
CHECKPOINT_FILE = "collatz_checkpoint.json"
RECORDS_FILE = "collatz_records.log"
ANOMALY_FILE = "collatz_ANOMALY.log"  # 異常（ループ・発散の疑い）を検知したときのログ
NICE_LEVEL = 19                 # プロセス優先度（Unix系。19が最も低い=他を邪魔しない）
MAX_STEPS = 10_000_000          # これを超えても1に到達しない場合は「異常」として扱い、無限ループを防ぐ
CACHE_MAX_SIZE = 500_000        # 軌道の「合流」を利用した高速化キャッシュの最大サイズ
DIGITS = len(str(START_NUMBER))  # 新しくランダムな数を選ぶときに使う桁数
_env_runtime = os.environ.get("COLLATZ_MAX_RUNTIME_SECONDS")
MAX_RUNTIME_SECONDS = float(_env_runtime) if _env_runtime else None  # Noneなら無制限
KTABLE_BITS = 18                # kステップテーブルの k。大きいほど速いが構築コストとメモリが増える
# ================================================================================

_CACHE = {}   # value -> (そこから1に到達するまでのステップ数, その間の最大値)
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


def random_start():
    """DIGITS桁のランダムな奇数を1つ生成する"""
    n = random.randint(10 ** (DIGITS - 1), 10 ** DIGITS - 1)
    return n | 1  # 奇数にする


def handle_stop(signum, frame):
    global running
    print("\n停止シグナルを受け取りました。安全に終了します...")
    running = False


signal.signal(signal.SIGINT, handle_stop)
signal.signal(signal.SIGTERM, handle_stop)


def lower_priority():
    """OSレベルでプロセスの優先度を下げる（できる範囲で）"""
    try:
        os.nice(NICE_LEVEL)  # Unix / macOS / Linux
        print(f"プロセス優先度を下げました (nice={NICE_LEVEL})")
    except (AttributeError, OSError):
        # Windows など os.nice が使えない環境
        try:
            import psutil
            p = psutil.Process(os.getpid())
            if sys.platform == "win32":
                p.nice(psutil.IDLE_PRIORITY_CLASS)
                print("プロセス優先度を下げました (IDLE_PRIORITY_CLASS)")
        except ImportError:
            print("優先度調整はスキップしました（psutil未インストール）。"
                  "動作には影響ありません。")


def collatz_steps(n, max_steps=MAX_STEPS):
    """コラッツ数列が1に到達するまでのステップ数と、途中の最大値を返す。

    3つの高速化を組み合わせている:
    1. 合流キャッシュ: 近くの数どうしは軌道が「合流」しやすい性質を利用し、
       一度計算した軌道の後半部分をキャッシュして再利用する。
    2. kステップテーブル: n mod 2^k だけで、次のkステップ分の変換を
       1回の掛け算・足し算・シフトにまとめて一気に進める。
    3. gmpy2: 利用可能なら、巨大整数演算そのものを高速なライブラリに委譲する。

    max_stepsを超えても1に到達しない場合は、無限ループを避けるため打ち切り、
    (steps, max_val, hit_limit=True) を返す。通常時は hit_limit=False。
    """
    if n == 1:
        return 0, 1, False

    visited = []  # (value, そのジャンプで進んだステップ数) のリスト
    m = _mpz(n)
    steps = 0
    while True:
        if m == 1:
            rem_steps, rem_max = 0, 1
            break
        m_key = int(m)
        if m_key in _CACHE:
            rem_steps, rem_max = _CACHE[m_key]
            break

        if m >= _KMOD:
            r = int(m & (_KMOD - 1))
            mult, add = _KTABLE[r]
            m_next = (mult * m + add) >> KTABLE_BITS
            step_inc = KTABLE_BITS
        else:
            m_next = m // 2 if m % 2 == 0 else 3 * m + 1
            step_inc = 1

        visited.append((m_key, step_inc))
        m = m_next
        steps += step_inc
        if steps >= max_steps:
            return steps, max(v for v, _ in visited) if visited else int(m), True

    suffix_steps = rem_steps
    suffix_max = rem_max
    for v, inc in reversed(visited):
        suffix_steps += inc
        suffix_max = max(suffix_max, v)
        if len(_CACHE) < CACHE_MAX_SIZE:
            _CACHE[v] = (suffix_steps, suffix_max)

    return steps + rem_steps, suffix_max, False


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            data = json.load(f)
        print(f"チェックポイントを読み込みました: 最長記録={data['best_steps']}ステップ, "
              f"これまでに{data['total_checked']:,}個チェック済み")
        return data
    return {
        "best_steps": 0,
        "best_n": None,
        "total_checked": 0,
    }


def save_checkpoint(data):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(data, f)


def log_record(n, steps, max_val):
    with open(RECORDS_FILE, "a") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | "
            f"steps={steps} | digits={len(str(n))} | "
            f"max_val_digits={len(str(max_val))} | n={n}\n"
        )


def log_anomaly(n, steps):
    with open(ANOMALY_FILE, "a") as f:
        f.write(
            f"{time.strftime('%Y-%m-%d %H:%M:%S')} | "
            f"n={n} が {steps:,}ステップを超えても1に到達しませんでした。"
            f"ループまたは発散の可能性があります。要確認。\n"
        )


def main():
    lower_priority()
    state = load_checkpoint()
    best_steps = state["best_steps"]
    start_time = time.time()

    print("探索を開始します。停止するには Ctrl+C を押してください。")
    if MAX_RUNTIME_SECONDS:
        print(f"制限時間: {MAX_RUNTIME_SECONDS/3600:.1f}時間で自動的に安全終了します。")
    print(f"毎回{DIGITS}桁のランダムな奇数を選んで探索します。")

    global running
    while running:
        if MAX_RUNTIME_SECONDS and time.time() - start_time >= MAX_RUNTIME_SECONDS:
            print(f"制限時間({MAX_RUNTIME_SECONDS/3600:.1f}時間)に達したため、安全に終了します。")
            running = False
            break

        for _ in range(BATCH_SIZE):
            if not running:
                break

            n = random_start()
            steps, max_val, hit_limit = collatz_steps(n)

            if hit_limit:
                # MAX_STEPSを超えても1に到達しなかった → ループ/発散の疑い
                # 無限ループにせず、ここで安全に停止して調査できるようにする
                log_anomaly(n, steps)
                print("=" * 60)
                print(f"⚠️  異常検知: n={n} が {steps:,}ステップを超えても収束しません")
                print(f"    {ANOMALY_FILE} に記録し、安全のため探索を停止します。")
                print("=" * 60)
                save_checkpoint(state)
                return

            state["total_checked"] += 1

            if steps > best_steps:
                best_steps = steps
                state["best_steps"] = best_steps
                state["best_n"] = n
                log_record(n, steps, max_val)
                print(f"新記録！ n={n} ({len(str(n))}桁) -> {steps}ステップ")

        save_checkpoint(state)

        # ここで意図的に休憩を入れてCPU負荷を抑える
        time.sleep(SLEEP_SECONDS)

    save_checkpoint(state)
    print(f"終了しました。（チェックしたのは合計{state['total_checked']:,}個）")


if __name__ == "__main__":
    main()
