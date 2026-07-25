#!/usr/bin/env python3
"""
fit_gamma.py

Post-processes the CSV files produced by coupling_scaling.c
(one file per bit-length window: <prefix>_bits<L>.csv, columns
agree,merge_time) and fits an exponential decay rate gamma to the
survival function P(tau_couple > t | early-merged) within each
window, to see how gamma(bits) drifts with scale.

Usage:
    python3 fit_gamma.py <prefix> <bit1> <bit2> ...

Example:
    python3 fit_gamma.py run1 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536
"""
import sys
import csv
import math
from collections import Counter


def load_merge_times(path):
    times = []
    agree = 0
    total = 0
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            total += 1
            a = int(row["agree"])
            mt = int(row["merge_time"])
            if a == 1:
                agree += 1
                if mt >= 0:
                    times.append(mt)
    return times, agree, total


def fit_exponential_tail(times, t_lo=5, t_hi_frac=0.8):
    """
    Fit log(survival(t)) ~ -gamma * t + c over an intermediate range
    of t (skip the very start where discreteness dominates, and the
    very end where counts are too sparse).
    Returns (gamma, r2, n_points) or None if not enough data.
    """
    if len(times) < 30:
        return None
    c = Counter(times)
    max_t = max(times)
    n = len(times)
    cum = 0
    survival = []
    for t in range(0, max_t + 2):
        cum += c.get(t, 0)
        survival.append((n - cum) / n)

    t_hi = max(t_lo + 5, int(max_t * t_hi_frac))
    xs, ys = [], []
    for t in range(t_lo, min(t_hi, len(survival))):
        if survival[t] > 0:
            xs.append(t)
            ys.append(math.log(survival[t]))
    if len(xs) < 5:
        return None

    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    varx = sum((x - mx) ** 2 for x in xs)
    if varx == 0:
        return None
    slope = cov / varx
    intercept = my - slope * mx

    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return (-slope, r2, len(xs))


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    prefix = sys.argv[1]
    bit_windows = [int(b) for b in sys.argv[2:]]

    print(f"{'bits':>8} {'samples':>8} {'agree_frac':>11} {'early_merge_frac':>17} "
          f"{'mean_tau':>10} {'gamma':>10} {'R2':>8}")
    for bits in bit_windows:
        path = f"{prefix}_bits{bits}.csv"
        try:
            times, agree, total = load_merge_times(path)
        except FileNotFoundError:
            print(f"{bits:>8}  (missing: {path})")
            continue

        agree_frac = agree / total if total else float("nan")
        early_merge_frac = len(times) / agree if agree else float("nan")
        mean_tau = sum(times) / len(times) if times else float("nan")

        fit = fit_exponential_tail(times)
        if fit:
            gamma, r2, npts = fit
            print(f"{bits:>8} {total:>8} {agree_frac:>11.4f} {early_merge_frac:>17.4f} "
                  f"{mean_tau:>10.3f} {gamma:>10.5f} {r2:>8.4f}")
        else:
            print(f"{bits:>8} {total:>8} {agree_frac:>11.4f} {early_merge_frac:>17.4f} "
                  f"{mean_tau:>10.3f} {'--':>10} {'--':>8}")


if __name__ == "__main__":
    main()
