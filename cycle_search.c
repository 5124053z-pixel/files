"""
Main demo: repeatedly prepend a fixed block to a number and track how the
Collatz total stopping time changes, showing convergence of the step
difference to the block's own bit length.

Usage:
    python3 collatz_block_repeat.py [x] [num_iterations]

Example:
    python3 collatz_block_repeat.py 10251997 300
"""
import sys
from collections import Counter


def collatz_steps(n: int) -> int:
    """Collatz total stopping time: steps until n reaches 1."""
    steps = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        steps += 1
    return steps


def run(x: int, num_iters: int = 300):
    L = x.bit_length()
    N = x
    prev = collatz_steps(N)
    diffs = []
    for k in range(2, num_iters):
        N = (x << N.bit_length()) | N
        s = collatz_steps(N)
        diffs.append(s - prev)
        prev = s
    return L, diffs


def main():
    x = int(sys.argv[1]) if len(sys.argv) > 1 else 10251997
    num_iters = int(sys.argv[2]) if len(sys.argv) > 2 else 300

    L, diffs = run(x, num_iters)
    print(f"x = {x} (L = {L} bits)")
    print(f"Ran {len(diffs)} iterations.\n")

    window = 50
    print(f"{'window':>8} {'freq(diff==L)':>15} {'freq %':>8}")
    for i in range(0, len(diffs), window):
        chunk = diffs[i:i + window]
        if len(chunk) < window:
            continue
        cnt = sum(1 for d in chunk if d == L)
        print(f"{i // window + 1:>8} {cnt:>15} {100 * cnt / len(chunk):>7.1f}%")

    overall = sum(1 for d in diffs if d == L) / len(diffs) * 100
    print(f"\nOverall freq(diff == L={L}): {overall:.1f}%")

    c = Counter(diffs)
    print(f"Top diff values: {c.most_common(5)}")


if __name__ == "__main__":
    main()
