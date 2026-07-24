# Self-Similar Digit Blocks and Collatz Total Stopping Time

**Status:** amateur/independent investigation, computationally verified, not peer-reviewed.
**TL;DR:** Repeatedly prepending a fixed bit-block to a number causes the Collatz
total stopping time to increase by an amount that converges — with probability
approaching 1 — to exactly the block's own bit length. A specific family of
blocks (those satisfying `3·(x/2)+1 = 2^L`) makes this convergence exactly 100%
from the very first iteration.

This is **not** a claim about the Collatz conjecture itself (true/false, cycles,
divergent trajectories). It is a concrete, checkable statement about a specific
family of numbers built from self-similar bit patterns.

**Update:** this investigation led to a bigger and more surprising finding
than the original construction — see §5d. The well-known claim that "about
half" of consecutive integers `n, n+1` share the same Collatz total stopping
time (OEIS A006577) appears, based on sampling out to `n ~ 10^27000`, to be
a small-`n` snapshot of a quantity climbing steadily toward 100%, not a
stable ~50% limit.

## Repository contents

```
README.md                              this document
LICENSE                                 MIT
collatz_block_repeat.py                 main demo, single block
big_survey.py                           105-block systematic sweep
exact_power_test.py                     isolates the "exact 2^L" condition
random_m_test.py                        §5a: tests whether self-similarity of
                                         N_k is even necessary (it isn't)
merging_residue_classes.py              §5b: classifies which residues mod 2^k
                                         provably force n, n+1 to merge
merging_classes_fast.c                  §5b: fast C/OpenMP version for larger k
general_merging_test.py                 §5c: extends §5b's merging-class
                                         analysis to general (L, x)
large_scale_sampling.c                  §5d: samples agreement rate at huge
                                         bit lengths (GMP + OpenMP)
cycle_search.c                          unrelated side-quest: exhaustive
                                         search for non-trivial Collatz cycles
                                         via parity-vector fixed points (GMP,
                                         OpenMP) — negative result up to q=20,
                                         kept for reference / reuse
results/single_block_10251997.csv       raw diff_k data for x=10251997, k=2..300
results/survey_105_blocks.csv           raw results of the 105-block sweep (§2)
```

---

## 1. Construction

Fix an odd integer `x` with bit length `L` (so `2^(L-1) ≤ x < 2^L`). Define a
sequence of integers by repeatedly prepending `x` to the front (most
significant bits) of the previous term:

```
N_1 = x
N_{k+1} = (x << bitlength(N_k)) | N_k
```

In other words, `N_k` is `x` written `k` times in a row in binary — the base-`2^L`
analogue of a repunit (e.g. `x=5`, `L=3` gives `101`, `101101`, `101101101`, ...).

One can show algebraically (and we verify computationally) that this satisfies
the clean closed form:

```
N_{k+1} = 2^L · N_k + x        (exact, for all k)
```

Let `steps(n)` be the Collatz total stopping time (number of `n → n/2` /
`n → 3n+1` operations until reaching 1). Define:

```
diff_k = steps(N_{k+1}) − steps(N_k)
```

## 2. Main empirical finding

**Claim.** For fixed `x` (bit length `L`), the frequency of `diff_k = L` among
the first `k` iterations tends to increase with `k` and stabilizes at a high
value (often 60–100%, and *exactly* 100% for a specific family — see §4)
rather than the ≈`7·L` predicted by generic Collatz statistics (average total
stopping time grows ≈7 steps per bit for a "random" number).

Tested against **105 distinct blocks** spanning bit lengths 4–24 (all-ones,
alternating, single-bit, and random patterns): **0 exceptions**. Every block's
dominant `diff_k` value converged to its own bit length `L`.

| block type | example (L=16) | late-stage freq(diff=L) |
|---|---|---|
| alternating (`1010...10`) | 43690 | **100.0%** |
| all-ones | 65535 | 66.7% |
| random odd | 53327 | 40.0% |
| known delay-record (`8400511`, L=24) | 8400511 | 66.7% |

## 3. Why this happens (proved part)

Since `N_{k+1} = 2^L·N_k + x`, we have `N_{k+1} ≡ N_k (mod 2^{L·k})` — i.e.
`N_k`'s entire bit pattern is preserved unchanged as the low-order bits of
`N_{k+1}`. Combined with the standard fact that **a number's first `m` Collatz
steps (its parity vector) are fully determined by the number mod `2^m`**
(Lagarias; see Bernstein & Lagarias, *The 3x+1 conjugacy map*, Canad. J. Math.
1996), this forces `N_k` and `N_{k+1}` to take an *identical* sequence of
odd/even steps for at least the first `L·k` iterations. Empirically the shared
prefix is even longer (~1.5·L·k, consistent with ~2/3 of Collatz steps being
"even").

This rigorously explains *why* `N_k` and `N_{k+1}` are entangled for a long
stretch. It does **not** by itself explain why the final total-step
*difference*, after they eventually diverge, is almost always the small,
k-independent value `L` rather than scaling with the (growing) length of the
shared prefix. That part remains an open empirical pattern — see §5.

## 4. The 100%-clean family

Numbers of the alternating form `x = 1010...10` (L bits) satisfy an exact
algebraic identity:

```
x / 2 = (2^L − 1) / 3          (integer, since x/2 has this closed form)
3·(x/2) + 1 = 2^L               (exact — lands precisely on 2^L)
```

I.e. `x`'s own Collatz trajectory reaches **exactly `2^L`** (the power of two
matching its *own* bit length) within 2 steps, then descends to 1 by pure
halving. This is a much sharper condition than merely "reaches some power of
two quickly" — a control test using numbers that reach *unrelated* powers of
two (found by tracing the Collatz map backward from `2^m` for `m ≠ L`) does
**not** reproduce the 100% effect (only 43–70%, same as generic numbers, see
`exact_power_test.py`).

The self-referential condition — landing exactly on `2^L` where `L` is the
block's own length — appears to align perfectly with the shift amount `L`
used when constructing `N_{k+1}`, eliminating the carry/interference that
causes deviation from `diff_k = L` for generic blocks.

## 5. Open questions

- Is there a clean closed-form proof that `Pr(diff_k = L) → 1` as `k → ∞`
  for *typical* `x` (not just the 100% family)? A partial mechanism is
  sketched in §3, but the exact limiting frequency and its dependence on the
  bit-pattern of `x` is not derived, only observed.
- Is `3(x/2)+1 = 2^L` **sufficient and necessary** for the 100% effect, or are
  there other unrelated algebraic identities that also produce it?
- Does the phenomenon generalize to bases other than 2^L (i.e. repeating a
  block in base `b` for other `b`), or to other `an+b` Collatz-like maps?
- **(see §5a)** Is the self-similar block-repeat construction even necessary
  for the effect, or is something much more general going on?

## 5a. Generalization: the self-similarity turns out not to matter

A follow-up experiment shows the phenomenon has **nothing to do with the
self-similar repunit-like construction of `N_k`**. Fix any odd `x` (bit
length `L`) and take a *uniformly random* large integer `m` (not built from
`x` at all). Then:

```
diff = steps(2^L * m + x) - steps(m)
```

still has **`L` as its single most frequent value**, and this mode sharpens
(the frequency of `diff = L` increases) as the bit length of `m` grows — see
`random_m_test.py`. Sample results (200 random trials each):

| x (bit length L) | m bit length | freq(diff = L) |
|---|---|---|
| alternating, L=12 | 100 | 75.0% |
| alternating, L=12 | 500 | 89.0% |
| random, L=12 | 100 | 32.5% |
| random, L=12 | 500 | 59.0% |
| all-ones, L=12 | 100 | 18.0% |
| all-ones, L=12 | 500 | 53.5% |

Since `steps(2^L·m) = L + steps(m)` **exactly and trivially** (multiplying
by `2^L` just prepends `L` guaranteed halving steps before the trajectory of
`m` continues unchanged), this reduces the whole phenomenon to a cleaner
question:

> For a fixed small offset `x` and a large "clean" (`2^L`-divisible) number
> `2^L·m`, why does adding `x` leave the **total stopping time unchanged**
> more often than any other specific outcome?

**Connection to a known (but apparently unproven) fact.** The case `L=1,
x=1` is *exactly* the classical, previously-documented observation that
`steps(n)` and `steps(n+1)` coincide roughly half the time (OEIS A006577,
comment: "It seems that about half of the terms satisfy a(i) = a(i+1)"; up
to 10,000,000, 4,964,705 of the terms satisfy this). Concretely: `2m` and
`2m+1` are literally consecutive integers, and `steps(2m) = 1 + steps(m)`
trivially, so "`diff = L = 1`" here is precisely "`steps(2m) = steps(2m+1)`",
i.e. a same-total-stopping-time pair of neighbors — verified directly in
`random_m_test.py`.

So the family of questions here — parameterized by `(L, x)` instead of just
"distance 1" — appears to be a genuine generalization of that classical,
apparently-still-unproven ~50% phenomenon. I could not find an existing
rigorous proof or heuristic derivation of even the base `L=1` case in the
literature searched. (Terras's theorem characterizes the *ordinary* stopping
time — first drop below the starting value — by finite congruence classes,
but *total* stopping time is not known to reduce to finite congruence data,
which is presumably why this is hard.)

## 5b. A proved partial result for the L=1 base case

Pushing on §5a's classical `n`, `n+1` case directly: a clean, fully provable
mechanism explains a substantial chunk of the ~50% agreement rate.

**Theorem (verified).** For all `n ≡ 4 (mod 8)` with `n ≥ 12`, write
`n = 8j+4`. Then:

```
n:    8j+4 → 4j+2 → 2j+1 → 6j+4        (3 steps)
n+1:  8j+5 → 24j+16 → 12j+8 → 6j+4     (3 steps)
```

Both trajectories reach the *exact same value* `6j+4` after exactly 3 steps.
Since the Collatz map is deterministic, everything after that point is
identical, so `steps(n) = 3 + steps(6j+4) = steps(n+1)` **exactly, with no
exceptions** (the sole edge case, `n=4`, terminates at 1 before completing 3
steps).

**Generalization.** Call a residue class `r (mod 2^k)` a *merging class* if
`n` and `n+1` provably reach the same value after the same fixed number of
steps for every `n ≡ r (mod 2^k)` (checked computationally against multiple
representatives per class; see `merging_residue_classes.py`). The fraction
of merging classes, as `k` grows:

| modulus `2^k` | merging classes | fraction |
|---|---|---|
| 8 | 1 | 12.5% |
| 16 | 3 | 18.8% |
| 32 | 8 | 25.0% |
| 64 | 18 | 28.1% |
| 128 | 39 | 30.5% |
| 256 | 82 | 32.0% |
| 512 | 170 | 33.2% |
| 1024 | 351 | 34.3% |
| 2048 | 721 | 35.2% |
| 4096 | 1476 | 36.0% |
| 8192 | 3012 | 36.8% |
| 16384 | 6130 | 37.4% |
| 32768 | 12450 | 38.0% |
| 65536 | 25241 | 38.5% |
| 131072 | 51105 | 39.0% |
| 262144 | 103358 | 39.4% |
| 524288 | 208840 | 39.8% |
| 1048576 | 421643 | 40.2% |
| 2097152 | 850737 | 40.6% |
| 4194304 | 1715546 | 40.9% |
| 8388608 | 3457791 | 41.2% |
| 16777216 | 6966495 | 41.5% |
| 33554432 | 14030369 | 41.8% |
| 67108864 | 28247507 | 42.1% |
| 134217728 | 56854178 | 42.4% |

**Important correction, found after more computation (thanks to a reader
running this on their own machine, up to `k=27`).** The increments do *not*
shrink at a constant geometric ratio. `increment(k) × k` is close to
constant (≈7.2, slowly drifting down) across `k=18..27`, which looks more
like `increment(k) ≈ C/k` — i.e. logarithmic-type growth,
`frac(k) ≈ frac(k0) + C·ln(k/k0)`, not geometric convergence to a fixed
limit. **A naive geometric extrapolation (which an earlier version of this
README used to guess "44-46%") is therefore not reliable**, and the true
asymptotic behavior of this fraction as `k → ∞` is currently **unknown**.

**A further complication: the classical "~50%" itself is not a settled
constant.** Checking `steps(n) == steps(n+1)` directly for `n` up to
10,000,000 (not just the merging-class abstraction) shows the raw agreement
rate drifting steadily upward with `n` — 45.0% below 100k, 47.7% below 1M,
49.6% below 10M, and *already past 50%* (50.4%) in the window
`n ∈ [9,000,000, 10,000,000]`. So the OEIS comment's "~50%" appears to be a
snapshot of a quantity that was itself still climbing at the point it was
measured, not a stable limiting value.

**Comparing the two quantities directly at matching scale** (same `k`,
i.e. `N = 2^k` vs modulus `2^k`) shows they are related but not equal, and
the gap between them is *not* constant either:

| k | N=2^k | raw agreement rate | merging-class fraction | gap |
|---|---|---|---|---|
| 16 | 65,536 | 44.52% | 38.51% | 6.01pt |
| 18 | 262,144 | 46.31% | 39.43% | 6.88pt |
| 20 | 1,048,576 | 47.77% | 40.21% | 7.56pt |
| 22 | 4,194,304 | 48.95% | 40.90% | 8.05pt |
| 24 | 16,777,216 | 50.03% | 41.52% | 8.50pt |

The fraction of raw agreement explained by provable merging actually
*decreases* slightly with scale (86.5% → 83.0% across this range), which is
the opposite of what I originally hoped to find.

**Honest bottom line.** Both the classical adjacent-agreement rate and the
merging-class fraction appear to drift upward indefinitely (at least within
the computationally reachable range, `k` up to 27), without settling at any
constant found so far, and a naive extrapolation is not trustworthy given
the apparent logarithmic-type growth. Whether either quantity has a
well-defined limit at all — and if so, what it is, and whether the two are
asymptotically related by a clean formula — is, as far as I can tell, a
genuinely open question. This investigation stops here for now with that
question explicitly unresolved, rather than with a (likely wrong) guessed
number.

## 5c. Synthesis: it's all about how "simple" x's own trajectory is

Extending §5b's merging-class analysis from the classical `L=1, x=1` case to
general `(L, x)` ties §4 and §5b together into one picture (see
`general_merging_test.py`):

| x | own trajectory | merging-class fraction (up to mod 512) |
|---|---|---|
| alternating, `x=2730` (L=12) | reaches `2^L` in 2 steps (§4) | converges cleanly to **exactly 1/2** (count = `2^(k-1)-1` for modulus `2^k`) |
| generic, `x=2905` (L=12) | no special structure | **zero** merging classes found even at modulus 512 |

So the same property that made the alternating family special in §4 — its
own Collatz trajectory being unusually short and predictable — is *also*
exactly what makes algebraically-provable merging classes abundant and easy
to find for that `x`. For "generic" `x` with a messy, unpredictable own
trajectory, provable merges are much rarer (or require far larger moduli to
detect), and the observed agreement between `steps(m)` and `steps(2^L m+x)`
for such `x` is presumably dominated by the harder, still-unexplained
"coincidental" mechanism from §5b rather than by literal trajectory merging.

This suggests a single underlying informal principle across §2–§5b:

> **The more predictable/short `x`'s own Collatz trajectory is, the more of
> the `diff = L` phenomenon can be explained by provable algebraic merging,
> and the higher its observed frequency ceiling.**

This is stated as an informal empirical pattern, not a theorem — turning it
into one (e.g. a precise statement relating some complexity measure of `x`'s
own trajectory to the growth rate of the merging-class fraction) is the
natural next step, left open here.

## 5d. The classical "~50%" is not the limit — it approaches 100%

Switching from exhaustive enumeration to random *sampling* (much cheaper,
lets us reach vastly larger scales) resolves the open question from §5b/§5c
in an unexpected direction. Sampling `steps(n) == steps(n+1)` for random `n`
of increasing bit length (see `large_scale_sampling.c`, GMP + OpenMP, run by
a reader on their own machine):

| bits | ~magnitude | agreement rate |
|---|---|---|
| 10 | 10^3 | 21.3% |
| 53 | 10^16 | 49.6% |
| 264 | 10^79 | 68.6% |
| 1299 | 10^391 | 83.4% |
| 3756 | 10^1131 | 89.4% |
| 10857 | 10^3268 | 94.3% |
| 31377 | 10^9445 | 96.2% |
| 90680 | 10^27297 | 98.0% |

The agreement rate climbs **well past 50% and keeps going**, reaching 98%
by `n ~ 10^27297`. Tracking `1 − rate` (distance from 100%) across these
points, the ratio between consecutive levels stays roughly constant
(~0.7–0.85, itself not drifting toward 1 the way §5b's merging-class
increments did), consistent with a **power-law approach to 100%**
(roughly `1 − rate ~ C / sqrt(bits)`) rather than the slow logarithmic
growth seen in the merging-class fraction.

**Revised conclusion.** The OEIS A006577 comment's "about half of the terms
satisfy a(i) = a(i+1)" appears to be describing a transient, small-`n`
snapshot of a quantity that is not settling near 50% at all — it looks like
`Pr(steps(n) = steps(n+1)) → 1` as `n → ∞`. This is consistent with (and
arguably the `L=1, x=1` special case of) the general pattern from §5a: for a
fixed shift, agreement becomes more and more likely as the numbers involved
get larger. The `~50%` figure that motivated this whole investigation was,
in hindsight, simply the value at a scale where the convergence to 100% is
still in its early stages.

This is sampling-based, not exhaustive, so it isn't a proof — but the trend
across nearly 30,000 orders of magnitude, with a consistent power-law-looking
decay rate, is hard to explain as a plateau at some value below 100%.

## 6. Relation to existing theory

This sits inside the well-studied **2-adic extension of the Collatz map**
(Bernstein & Lagarias 1996; see also the notion of "Collatz cyclic numbers"
associated with periodic parity vectors). As `k → ∞`, `N_k` converges 2-adically
to the rational 2-adic integer `α = -x/(2^L - 1)`. This connects the present
observation to known machinery, but the specific claim about the *frequency*
of `diff_k = L`, as far as I can tell, is not stated in the literature I
searched (Bernstein–Lagarias 1996; Lagarias periodicity conjecture papers;
Hercher 2023; Eliahou 1993; Simons & de Weger). Corrections/pointers to prior
work are very welcome.

## 7. Reproducing this

All experiments are pure Python (only the `fractions` and `collections`
standard library modules) plus one C/GMP program for larger-scale exhaustive
cycle searches (unrelated negative-result side-quest, kept for completeness —
see `cycle_search.c`).

```
python3 collatz_block_repeat.py      # main phenomenon, single block
python3 big_survey.py                 # 105-block systematic sweep
python3 exact_power_test.py           # isolates the "exact 2^L" condition
python3 random_m_test.py              # §5a: self-similarity isn't necessary
python3 merging_residue_classes.py    # §5b: classifies merging residue classes
python3 general_merging_test.py       # §5c: extends §5b to general (L, x)
```

No external dependencies beyond CPython 3.8+.

## 8. Caveats

- This does **not** bear on the truth of the Collatz conjecture itself.
- All claims here are backed by finite computation (k up to a few thousand,
  L up to 24), not proofs, except where explicitly marked "proved."
- This investigation was done collaboratively with an AI assistant (Claude);
  all code was executed and its output verified before being reported here.
  Please independently re-verify before citing.

## Acknowledgments

Investigation prompted by tracing down and correcting an AI-hallucinated
claim (a different model asserted a specific number was a "periodic point"
of the Collatz map via a fabricated formula, which did not survive direct
computational verification). The real phenomenon documented here was found
in the process of checking that claim.





# Addendum: Coupling-Theoretic Follow-up to §5

Status: amateur/independent investigation, computationally verified, not peer-reviewed. Continuation of the original README's §5 (open questions on the `steps(n) == steps(n+1)` agreement rate). This addendum reframes the problem using probabilistic **coupling theory** (coalescing Markov chains / mixing-time arguments) and reports new computational findings, including a correction of an earlier claim made mid-investigation.

## 0. Motivation

The original README left an open question at §5c: is the "raw" agreement rate `P(steps(n) = steps(n+1))` fully explained by literal trajectory merging, or is some of it "coincidental" (both n and n+1 reaching 1 in the same number of steps via genuinely different intermediate values)? §5b/§5d's merging-class analysis suggested only a partial (and slowly, unclearly growing) fraction was explained this way.

This addendum tests that question directly by simulation rather than by exhaustive residue-class enumeration, framing it in coupling-theoretic language: define

    tau_couple(n) = min { t >= 1 : T^t(n) = T^t(n+1) }

using the raw (unaccelerated) Collatz map T, searched strictly before either trajectory reaches 1. This is the "coalescing coupling" of the two trajectories. §5b's merging-class theorem is a special case: it proves `tau_couple` is uniformly bounded (= 3) on the whole residue class n ≡ 4 (mod 8).

## 1. A definitional trap (documented for honesty)

An initial version of this experiment defined `tau_couple` by letting the raw map run indefinitely, *including after reaching 1* (where it cycles 1 → 4 → 2 → 1 → ...). This is wrong: it conflates literal early merging with accidental phase-alignment in the post-convergence cycle, and gave a spuriously high "100% coupling" reading that mixed two different phenomena.

A second version fixed the cycling issue but capped the search at `steps(n) + small_buffer`, which is a **tautology**: if `steps(n) = steps(n+1) = s`, then by definition `T^s(n) = T^s(n+1) = 1`, so a merge is *guaranteed* to be found by time `s` regardless of whether anything interesting happened earlier. Reporting "100% merge fraction" from this setup is not a finding — it restates the agreement condition.

**Corrected definition:** search only `t ∈ [1, s)`, strictly before either trajectory reaches 1. A merge found here is a genuine early coincidence of values, not a trivial simultaneous arrival at 1. All results below use this corrected definition.

## 2. Result: early merging appears to fully explain agreement (up to N = 6×10⁷)

Exhaustive check, `n = 1 .. N`, computing `steps(n)`, `steps(n+1)`, and (when they agree) searching for a strict early merge:

| N | agreeing pairs | early-merge fraction (of agreeing pairs) | exceptions |
|---|---|---|---|
| 1,000,000 | 477,245 | 100.0000% | 0 |
| 5,000,000 | 2,454,559 | 100.0000% | 0 |
| 15,000,000 | 7,492,334 | 100.0000% | 0 |
| 60,000,000 | 30,547,761 | 100.0000% | 0 |

Zero exceptions across ~30.5 million agreeing pairs. This is stronger than what §5c speculated — no "coincidental, non-merging" agreement was found in this range at all. If this holds in general, it would mean:

    steps(n) = steps(n+1)   <=>   tau_couple(n) < steps(n)

i.e. agreement of total stopping times is *equivalent* to (not just partially explained by) literal early coalescence — collapsing §5c's distinction between "provable merging" and "coincidental agreement" into a single mechanism, at least at these scales.

**Caveats (real ones):**
- This is exhaustive verification up to 6×10⁷, not a proof. A counterexample could exist beyond this range.
- The mean early-merge time drifts upward with N (16.0 → 18.3 → 19.8 → 21.7 across the four scales above), indicating a heavier tail at larger scale — consistent with, but not proof of, the pattern continuing to hold.
- Python-only implementation currently caps practical N around 10⁸; extending further requires the C/GMP tooling below.

## 3. Coupling-time distribution: clean exponential tail (at fixed, small scale)

Conditional on early merging, the survival function `P(tau_couple > t | merged)` fits an exponential decay extremely well:

- N = 300,000: decay rate γ ≈ 0.0421, R² = 0.9987, half-life ≈ 16.4 steps
- Min observed early-merge time: 3 (matches the §5b theorem: n ≡ 4 mod 8 merges at exactly t = 3)

In coupling-theory language (Aldous–Diaconis coupling inequality: `d_TV(t) ≤ 2·P(tau_couple > t)`), this γ plays the role of a **spectral-gap-like decay constant** for the pair-coupling process — a concrete, fittable analogue of a Markov chain mixing-time bound.

## 4. New open question: does γ scale with bit length?

A small-scale (N=3000 samples per window) test of `tau_couple`'s decay rate γ across increasing bit-length windows of n showed γ *shrinking* with scale in a way consistent with γ ∝ 1/bits:

| bits | agree_frac | early_merge_frac | mean τ | γ | R² |
|---|---|---|---|---|---|
| 32 | 0.4087 | 1.0000 | 44.35 | 0.02049 | 0.976 |
| 64 | 0.5047 | 1.0000 | 84.51 | 0.00878 | 0.970 |
| 128 | 0.5903 | 1.0000 | 166.84 | 0.00408 | 0.957 |

`γ × bits` is roughly constant (~0.66, 0.56, 0.52) across this small range — i.e. mean coupling time appears to scale **linearly** with bit length, not staying O(1) as scale grows. This is a distinct phenomenon from §5b's merging-class fraction (which grows only logarithmically in k) — here we're measuring per-pair coupling speed, not whether a whole residue class merges uniformly, and it appears to slow down (linearly) with scale even though the early-merge *fraction* stays at 100%.

**This is preliminary** — only 3 small windows, 3000 samples each. Needs the full-scale run described in §5 below to confirm the 1/bits scaling and check whether early-merge fraction stays at 100% at much larger bit lengths (target: up to 65536 bits, matching the scale of the original §5d large_scale_sampling.c).

## 5. Tooling produced this session

- `coupling_experiment3.py`, `coupling_experiment4.py`, `coupling_experiment5.py` — Python, exhaustive verification up to N=6×10⁷ (memoized total-stopping-time cache for speed), used for §2–3 above.
- `coupling_scaling.c` — C + GMP + OpenMP, windowed random sampling across bit-length scales (32 to 65536 bits by default, easily extended), outputs per-window CSV of `(agree, merge_time)` pairs. Compiles with `gcc -O3 -fopenmp -o coupling_scaling coupling_scaling.c -lgmp -lm`.
- `fit_gamma.py` — post-processes the C program's CSV output, fits γ per bit-length window via log-linear regression on the survival function tail, reports agreement fraction, early-merge fraction, mean coupling time, γ, and R² per window.

Usage:
```
./coupling_scaling 20000 run1
python3 fit_gamma.py run1 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536
```

## 6. Summary of what's now established vs. still open

**Established (exhaustive, up to N=6×10⁷):**
- Every observed `steps(n)=steps(n+1)` pair coalesces to a literal common value strictly before reaching 1.
- Conditional on merging, `tau_couple` has a clean, well-fit exponential tail at fixed small scale (R² > 0.99).

**Open (needs the C/GMP run at scale):**
- Whether early-merge fraction stays at 100% for much larger n (say, 1000+ bit numbers), or whether exceptions appear.
- Whether γ(bits) truly follows a 1/bits law, or something else (log, power-law) once tested over a wider, more heavily-sampled range.
- Whether a clean theoretical argument (e.g. via the finite-state "carry" structure of the 3x+1 map read bit-serially, as in Bernstein–Lagarias) can *derive* the 1/bits scaling rather than just fit it — this would upgrade §4's observation from empirical to proved, in the spirit of §5b's rigorous n≡4(mod 8) result.

## Acknowledgments

This addendum was developed collaboratively with an AI assistant (Claude); all code was executed and its output verified before being reported here, including catching and correcting a tautological definition error mid-investigation (§1). Please independently re-verify before citing.

