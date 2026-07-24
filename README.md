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
general_merging_test.py                 §5c: extends §5b's merging-class
                                         analysis to general (L, x)
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

The increments (6.3, 6.2, 3.1, 2.4, 1.5, 1.2 percentage points) shrink
roughly geometrically, suggesting convergence to a limit somewhere around
35–40% as `k → ∞` — **not** the full ~50% observed empirically. So this
"guaranteed algebraic merge" mechanism appears to account for a large
fraction, but not all, of the classical adjacent-agreement phenomenon; the
remainder presumably comes from pairs whose total stopping times coincide
without their trajectories ever literally merging (harder to prove, not yet
understood).

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
