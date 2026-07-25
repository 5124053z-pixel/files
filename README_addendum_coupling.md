# Addendum: Coupling-Theoretic Follow-up to §5

Status: amateur/independent investigation, computationally verified, not peer-reviewed. Continuation of the original README's §5 (open questions on the `steps(n) == steps(n+1)` agreement rate). This addendum reframes the problem using probabilistic **coupling theory** (coalescing Markov chains / mixing-time arguments) and reports new computational findings, including corrections of errors made mid-investigation.

## 0. Motivation

The original README left an open question at §5c: is the "raw" agreement rate `P(steps(n) = steps(n+1))` fully explained by literal trajectory merging, or is some of it "coincidental" (both n and n+1 reaching 1 in the same number of steps via genuinely different intermediate values)? §5b/§5d's merging-class analysis suggested only a partial (and slowly, unclearly growing) fraction was explained this way.

This addendum tests that question directly by simulation rather than by exhaustive residue-class enumeration, framing it in coupling-theoretic language: define

    tau_couple(n) = min { t >= 1 : T^t(n) = T^t(n+1) }

using the raw (unaccelerated) Collatz map T, searched strictly before either trajectory reaches 1. This is the "coalescing coupling" of the two trajectories. §5b's merging-class theorem is a special case: it proves `tau_couple` is uniformly bounded (= 3) on the whole residue class n ≡ 4 (mod 8).

## 1. A definitional trap (documented for honesty)

An initial version of this experiment defined `tau_couple` by letting the raw map run indefinitely, *including after reaching 1* (where it cycles 1 → 4 → 2 → 1 → ...). This is wrong: it conflates literal early merging with accidental phase-alignment in the post-convergence cycle, and gave a spuriously high "100% coupling" reading that mixed two different phenomena.

A second version fixed the cycling issue but capped the search at `steps(n) + small_buffer`, which is a **tautology**: if `steps(n) = steps(n+1) = s`, then by definition `T^s(n) = T^s(n+1) = 1`, so a merge is *guaranteed* to be found by time `s` regardless of whether anything interesting happened earlier. Reporting "100% merge fraction" from this setup is not a finding — it restates the agreement condition.

**Corrected definition:** search only `t ∈ [1, s)`, strictly before either trajectory reaches 1. A merge found here is a genuine early coincidence of values, not a trivial simultaneous arrival at 1. All results below use this corrected definition.

## 2. Result: early merging appears to fully explain agreement, at every scale tested

**Exhaustive check** (Python), `n = 1 .. N`, computing `steps(n)`, `steps(n+1)`, and (when they agree) searching for a strict early merge:

| N | agreeing pairs | early-merge fraction | exceptions |
|---|---|---|---|
| 1,000,000 | 477,245 | 100.0000% | 0 |
| 5,000,000 | 2,454,559 | 100.0000% | 0 |
| 15,000,000 | 7,492,334 | 100.0000% | 0 |
| 60,000,000 | 30,547,761 | 100.0000% | 0 |

**Windowed random sampling** (C + GMP, `coupling_scaling.c`), single n per bit-length window:

| bits | samples | agree_frac | early_merge_frac | mean τ |
|---|---|---|---|---|
| 32 | 3000 | 40.87% | 100.0000% | 44.4 |
| 4096 | 2000 | 90.20% | 100.0000% | 2296.5 |
| 16384 | 800 | 95.63% | 100.0000% | 5011.9 |
| 65536 | 150 | 98.00% | 100.0000% | 13899.4 |

Zero exceptions across ~30.5 million exhaustively-checked pairs *and* every sampled bit-length window up to 65536 bits (n ~ 10^19728). No "coincidental, non-merging" agreement was found anywhere. The agreement-rate values also cross-check cleanly against the original README's §5d table (e.g. 98.0% at ~65536 bits, matching). If this holds in general:

    steps(n) = steps(n+1)   <=>   tau_couple(n) < steps(n)

i.e. agreement of total stopping times appears *equivalent* to (not just partially explained by) literal early coalescence — collapsing §5c's distinction between "provable merging" and "coincidental agreement" into a single mechanism, at every scale tested so far.

**Caveats:** this is exhaustive/sampled verification, not a proof; a counterexample could exist beyond the tested range. Mean early-merge time grows with scale (see §4 below), so the tail is getting heavier, not disappearing.

## 3. Coupling-time distribution: exponential-looking tail, but it's a two-component mixture

Conditional on early merging, a naive single-exponential fit to `P(tau_couple > t | merged)` looks clean at fixed small scale (N=300,000: γ≈0.0421, R²=0.9987), with minimum observed merge time 3 (matching the §5b theorem exactly). But this single-exponential picture turns out to be too simple — see §5.

## 4. γ(bits) scaling: roughly γ ∝ 1/bits, but with a caveat

Fitting the tail decay rate γ per bit-length window (log-linear regression on the survival function, windows 32–65536 bits):

| bits | γ | γ×bits |
|---|---|---|
| 32 | 0.02049 | 0.656 |
| 512 | 0.00098 | 0.502 |
| 4096 | 0.00013 | 0.532 |
| 16384 | 0.00003 | 0.492 |
| 65536 | 0.00001 | 0.655 |

Log-log regression of γ against bits gives exponent **β = −1.018** (the γ×bits column above is roughly flat across 3 orders of magnitude), consistent with γ ~ 1/bits.

However, fitting `mean(tau_couple)` against bits directly gives exponent **α = 0.728** (R² = 0.992) — not the α ≈ 1 a pure single-rate exponential model (mean = 1/γ) would predict. This mismatch was the motivation for §5.

## 5. Resolving the mismatch: the distribution is a two-component mixture

Quantile analysis (p10, p50, p90, p99, max of tau_couple) per bit-length window, each fit to a power law vs. bits:

| quantile | power-law exponent α | R² |
|---|---|---|
| p10 | **0.058** (essentially flat) | 0.738 |
| p50 (median) | 0.384 | 0.909 |
| p90 | 0.661 | 0.962 |
| p99 | **0.935** (nearly linear) | 0.998 |
| max | 0.951 | 0.9998 |

The exponent increases smoothly from ~0 to ~1 across the quantile range. Concretely, p10 is **numerically almost constant** (6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.8, 9.0) across the *entire* 32-to-65536-bit range — an 11-order-of-magnitude change in n's scale with essentially no change in the fastest 10%'s merge time. Meanwhile p99 and the max grow almost linearly with bits.

**Conclusion:** `tau_couple` is not a single exponential. It is a mixture of (a) a bits-independent fast component (roughly the bottom ~90% of the distribution, driven by small, fixed residue classes — see §6) and (b) a heavy tail whose scale grows roughly linearly with bit length. The previously-observed mean exponent (0.73) is what you get when averaging over a mixture whose bulk is O(1) and whose tail is O(bits); it sits strictly between the two component behaviors, as expected for a mixture, and is not itself a fundamental exponent.

## 6. Hypothesis 1 (confirmed): the fast component is explained by tiny, scale-independent residue classes

Direct test: for 150 "fast" mergers (tau_couple < 20) sampled at a fixed, large bit length (4096 bits), binary-search the smallest modulus 2^K such that fixing n's low K bits (with fresh random high bits) reliably reproduces the same merge time across multiple independent trials.

Result: **all 150/150 samples were fully explained by K ≤ 14** (mean K = 5.23):

| K | count |
|---|---|
| 3 | 48 |
| 4 | 24 |
| 5 | 29 |
| 6 | 15 |
| 7 | 11 |
| 8–14 | 23 |

K=3 (modulus 2^3=8) is the single most common value, at 48/150 (32%) — this is *exactly* the §5b theorem's residue class n ≡ 4 (mod 8), now shown empirically to be the dominant single mechanism behind the fast component of the distribution, not merely an isolated example. This directly confirms that §5's "fast, bits-independent" component is driven by small, fixed residue classes generalizing the §5b mechanism, matching the near-flat p10 behavior found in §5.

A related control test (fixing n's low 24 bits and checking whether merge time is reproduced) shows a sharp transition: merges with tau_couple below ~24–40 are ~100% reproduced by 24 bits of information; merges above that are essentially never reproduced. This is expected given the classical fact that a trajectory's first t steps are determined by n mod 2^t (Bernstein–Lagarias) — merges taking t steps generically need on the order of t bits to pin down — but confirms the mechanism is tight (no large hidden dependence on far-away bits for the fast component).

## 7. Hypothesis 2 (mostly refuted, with a clean positive finding instead): n's *own* trajectory length is irrelevant; local 2-adic structure is what matters

Tested whether §5c's principle ("the simpler x's own trajectory, the more merging is explained") transfers to individual n — i.e. does a short/simple own total stopping time `steps(n)` predict fast merging with n+1? Sampled 4000 pairs at 2048 bits, all with `steps(n)=steps(n+1)`:

- `corr(steps(n)/bits, log(tau_couple)) = +0.020` — **essentially zero**. Fastest/median/slowest 10% groups all have the same mean `steps(n)/bits ≈ 7.2`, matching the well-known generic ~7 steps/bit average with no distinguishing signal.
- Instead, **local 2-adic structure at the low-order bits** correlates meaningfully: `corr(v2(n+1), log τ) = +0.363`, `corr(run_length(n+1), log τ) = +0.289` (v2 = 2-adic valuation / trailing-zero count).
- The fastest 10% of all samples had **exactly** v2(n) = 2.000 (zero variance) and merge_time exactly 3 — i.e. the fastest decile is entirely the n ≡ 4 (mod 8) class, consistent with §6.

**An apparent asymmetry** — v2(n+1) correlating much more strongly than v2(n) (+0.363 vs. −0.021 pooled over all samples) — turned out to be a **pooling artifact**, not a real effect: since exactly one of n, n+1 is even, v2(n)=0 whenever n is odd and vice versa, so pooling both cases together dilutes whichever variable is "trivially zero" in each half. Comparing like-for-like (within the n-even subgroup only, where v2(n) is meaningful, vs. within the n-odd subgroup only, where v2(n+1) is meaningful) gives symmetric correlations: **+0.322 and +0.306 respectively** — the same effect, as expected by the n ↔ n+1 symmetry of the setup.

Partial correlation of `steps(n)/bits` against `log(tau_couple)`, controlling for v2(n+1) (and separately for run-length), stays near zero (+0.016, +0.015) — confirming the null result is not merely masked by the v2 confound. Stratifying by v2(n+1) value and computing within-stratum correlation of `steps(n)` with merge time also stays near zero in every stratum.

**Conclusion:** merging speed is governed by *local* low-order bit structure (2-adic valuation / run-length near the LSB, i.e. exactly the kind of information the Lagarias congruence fact says determines early steps), not by any *global* property of n's own trajectory shape. §5c's "simple trajectory" principle, formulated for the specific repeated-block construction of the original README (§1–§4), does not appear to generalize to arbitrary (n, n+1) pairs in the form tested here.

## 8. Tooling produced this session

- `coupling_experiment3.py` / `4.py` / `5.py` — Python, exhaustive verification up to N=6×10⁷ (memoized total-stopping-time cache), used for §2.
- `coupling_scaling.c` — C + GMP + OpenMP, windowed random sampling across bit-length scales (32–65536 bits, extensible), outputs per-window CSV of `(agree, merge_time)`. Build: `gcc -O3 -fopenmp -o coupling_scaling coupling_scaling.c -lgmp -lm`.
- `fit_gamma.py` — fits γ per bit-length window from the C program's output (log-linear regression on the survival-function tail).
- `fast_slow_test.py` — tests whether merges below/above a threshold are reproduced by fixing n's low K bits (§6, sharp-transition control test).
- `minimal_K_test.py` — binary search for the smallest modulus 2^K explaining each fast merger (§6 main result).
- `hypothesis2_test.py` — correlates merge time against n's own trajectory length and 2-adic/run-length features at fixed bit length; produces the partial-correlation and stratified-correlation analysis (§7).

## 9. Summary: established vs. still open

**Established (exhaustive up to N=6×10⁷; sampled up to 65536 bits with zero exceptions):**
- Every observed `steps(n)=steps(n+1)` pair coalesces to a literal common value strictly before reaching 1, at every scale tested.
- `tau_couple` is a two-component mixture: a bits-independent fast bulk (~90%, explained by tiny fixed residue classes, dominated by the §5b n≡4(mod 8) mechanism) plus a heavy tail growing roughly linearly in bit length (~top 1%, p99 exponent ≈ 0.94).
- Merging speed depends on *local* 2-adic structure near the LSB of n (and n+1), not on any global property of n's own Collatz trajectory length; the apparent v2(n) vs v2(n+1) asymmetry is a pooling artifact and disappears under like-for-like comparison.

**Open:**
- Whether the 100% early-merge / zero-exception pattern holds for all n (currently: exhaustive to 6×10⁷, sampled with zero exceptions to 65536 bits).
- A precise theoretical model for the two-component mixture — e.g. can the fast/slow split point and the p99 ≈ bits^0.94 exponent be derived (rather than fit) from the finite-state "carry" structure of the 3x+1 map read bit-serially (Bernstein–Lagarias)?
- Why the fast component's modulus needs (§6, K values 3–14) cluster where they do, and whether this connects to §4's "exact 2^L" family in the original README.

## Acknowledgments

This addendum was developed collaboratively with an AI assistant (Claude); all code was executed and its output verified before being reported here, including catching and correcting a tautological definition error mid-investigation (§1) and a data-construction bug that initially produced a misleading result (§6). Please independently re-verify before citing.
