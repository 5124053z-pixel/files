# Collatz search: what the data actually says

Analysis notes for the searches in this repository, written up after the
1997-chain log reached 26,060 iterations (104,240 digits).

Everything below is measured from `collatz_1997_records.log` or recomputed
from scratch. Where an earlier conclusion turned out to be wrong, the wrong
version and the reason it was wrong are kept, because the mistakes were more
instructive than the result.

## Summary

| Question | Answer |
|---|---|
| Does the 1997 chain show unusual behavior? | Yes, but not because it is 1997 |
| Is the effect real? | In decimal, weakly (z = −3.28 at best, and it vanishes if size is matched by bit length instead). In binary, strongly |
| Is it specific to 1997? | No. It appears in repeated-digit numbers generally |
| Do slow-falling families exist? | Yes, but only visible in base 2, not base 10 |
| Binary or decimal repetition more anomalous? | Binary, by roughly 10x in effect size |
| Is the cause understood? | **No** |
| Does it bear on the Collatz conjecture? | **No** — see "Why this is not a lead" |

## The observable

For a starting number `n`, let `steps(n)` be the number of Collatz steps to
reach 1, counting both `3n+1` and `n/2` as one step each, and let
`d = len(str(n))`. The quantity of interest is `steps(n) / d`.

The heuristic prediction is

```
steps/digit = 3 / log10(4/3) = 24.0118
```

which comes from: one odd step is followed on average by two divisions, and
each such cycle multiplies `n` by roughly `3/4`, so it takes `1/log10(4/3)`
cycles of 3 steps each to shed one decimal digit.

Uniformly random odd numbers match this. Measured means, with the sample
mean's own error bar:

| digits | samples | mean | 95% CI |
|---|---|---|---|
| 200 | 150 | 23.8472 | [23.56, 24.13] |
| 800 | 60 | 24.1352 | [23.90, 24.37] |
| 3000 | 20 | 24.0127 | [23.82, 24.21] |

## The finding

Numbers built by repeating a fixed digit block fall *faster* than random
numbers of the same size. At 600 digits, against a random baseline of
24.0585 (n = 40):

| block width | mean steps/digit | seeds | z |
|---|---|---|---|
| random | 24.0585 | 40 | — |
| 1 | 21.9174 | 9 | −6.39 |
| 2 | 23.7702 | 16 | −1.15 |
| 3 | 23.3195 | 16 | −2.94 |
| 4 | 23.0122 | 16 | −4.17 |
| 5 | 23.9376 | 16 | −0.48 |
| 6 | 23.7581 | 16 | −1.20 |

The effect is strong at widths 1, 3, and 4, and absent at 2, 5, and 6. It is
**not monotonic in width**, which rules out the simplest explanations.

Checked with sample sizes matched (8 repeated-digit numbers against repeated
draws of 8 random numbers, at 1200 digits):

```
repeated-digit, 8 seeds : 23.3536
random, 8 draws x 20    : 24.0017 +/- 0.1843
                      z : -3.52
```

The 1997 chain sits inside this family and is not distinguished within it.
Across widths and seeds it is unremarkable; at some sample points it is the
*highest* of the seeds tested, not the lowest.

This decimal effect turns out to be weaker than it looks here, and partly an
artifact of holding *decimal* length fixed — see "Which construction is more
anomalous" below, where matching on bit length instead drops it to z = −0.64.

### Mechanism: partially identified, not explained

The proximate cause is the ratio `r` = divisions per odd step, which feeds

```
steps/digit = (1 + r) / (r*log10(2) - log10(3))
```

This is steep: r = 2.000 gives 24.012, r = 2.010 gives 23.525, r = 2.020
gives 23.061. Measured, averaged over several sizes:

```
repeated-digit chains : r = 2.01286
random odd numbers    : r = 2.00199
```

A shift of +0.011 in r accounts for the observed gap. **Why repeated-digit
numbers have elevated r is not known.** Two candidate explanations were
tested and both failed:

- **Leading digits.** `len(str(n))` is a coarse stand-in for `log10(n)`, and
  repeated-digit numbers have a fixed `log10` fractional part while random
  ones follow Benford. Recomputing with `log10(n)` as the denominator moves
  the gap from −0.6218 to −0.6163 — essentially unchanged. Not the cause.
- **Residue mod 4.** Every term of the 1997 chain is ≡ 1 (mod 4), so `3n+1`
  is divisible by 4 and the first step always divides at least twice (the
  chain's first exponent is always exactly 3, against 1.959 for random odds).
  But this is one step out of several thousand, contributing ~1/6000 of r.
  Sampling random numbers restricted to n ≡ 1 and n ≡ 3 (mod 4) gives 24.0602
  and 24.1928 — both at the random baseline. Within a few steps the orbit's
  mod-4 distribution equalizes to ~0.5 either way. Not the cause.

### The symmetric question: can a structure fall *slower*?

Worth asking, since if elevated r makes numbers fall faster, suppressed r
should make them fall slower, and a family that reliably lags the baseline
would be the more interesting direction — slow-falling numbers are at least
adjacent to the idea of a number that never falls.

The attempt: numbers ≡ 3 (mod 4) have `3n+1` divisible by 2 exactly once,
which should depress r. Measured at 800 digits:

```
n = 3 (mod 4)  : 24.1928   r = 1.99709
n = 1 (mod 4)  : 24.0602   r = 1.99963
random         : 24.1357
```

Both sit at the baseline. Constraining the *starting* residue does not move
the ratio, for the same reason it failed to explain the repeated-digit
effect: the orbit's residues equalize within a few steps, so a condition on
`n` alone washes out.

No slow-falling family was found *in base 10*, and every significant
deviation among decimal block widths was negative. That turned out to be an
artifact of the base — see below.

## Binary is the right base, and it produces slow-falling families

Repeating a digit block in base 10 is the wrong construction: the map divides
by 2, so structure should be imposed in base 2. Redone with `n` built by
repeating a `w`-bit pattern to 1200 bits, and the ratio measured against
`log2(n)` rather than a decimal digit count:

```
steps/bit for a random orbit = 3 / log2(4/3) = 7.2283
```

Measured, with the degenerate all-ones pattern excluded (see caveat below),
against a random baseline of 7.1149 (n = 40, sigma = 0.2918):

| width w | mean | median | seeds | z | verdict |
|---|---|---|---|---|---|
| random | 7.1149 | 7.1030 | 40 | — | — |
| 3 | 7.8044 | 7.7973 | 6 | +5.79 | **slow** |
| 4 | 5.7978 | 5.3841 | 14 | −16.89 | **fast** |
| 5 | 7.7026 | 7.6663 | 20 | +9.01 | **slow** |
| 6 | 7.0835 | 7.4896 | 18 | −0.46 | none |
| 8 | 7.4333 | 7.8297 | 18 | +4.63 | **slow** |
| 10 | 8.2120 | 8.2252 | 20 | +16.82 | **slow** |
| 12 | 7.3378 | 7.2836 | 20 | +3.42 | **slow** |
| 16 | 7.0866 | 7.0519 | 17 | −0.40 | none |
| 20 | 7.1643 | 7.2246 | 16 | +0.68 | none |
| 24 | 7.1401 | 6.9595 | 18 | +0.37 | none |

Measured separately at 1206/1210 bits so the width divides evenly:
w = 9 gives 6.7617 (z = −4.35, **fast**) and w = 11 gives 7.1984
(z = −0.33, none).

So slow-falling families do exist — widths 3, 5, 8, 10 and 12 all sit
significantly above the random baseline, with w = 10 at z = +16.82. The
base-10 experiment found only fast families because repeating decimal digits
imposes no clean 2-adic structure; in binary both signs appear immediately.

Width 4 is the striking outlier in the other direction, and its median
(5.3841) is well below its mean, so it is not one stray value driving it.

**Caveat on two degenerate patterns.** The all-ones pattern gives
`2^1200 − 1` regardless of `w`, so it was being counted once per width as if
it were a different number each time, at 13.1517 — inflating every mean.
It is excluded above. Separately, `p = 0101` at w = 4 gives exactly
`(2^1200 − 1)/3`, ratio 1.0022; it is genuine, not an artifact, and is part
of why width 4 is extreme.

### The two outliers are explained, and they are the same phenomenon

Both extremes have exact closed-form reasons, and they sit at opposite ends
of one axis: how long the orbit can be forced to stay odd.

**`0101…01` = `(2^k − 1)/3`, the fastest possible.** For this n,
`3n + 1 = 2^k` exactly — verified as an integer identity at k = 1200. The
orbit therefore takes one odd step to land precisely on a power of two, then
halves k times and is done: `1 + 1200 = 1201` steps, matching the measured
value. Nothing about the trajectory is random; it is the shortest route a
1200-bit number can take.

**`111…1` = `2^k − 1`, the slowest.** Here `3n + 1 = 3·2^k − 2 =
2(3·2^{k−1} − 1)`, so exactly one division returns an odd number, and the
pattern repeats. Each odd-plus-one-division cycle multiplies n by about 3/2,
so the number *grows*. Traced over the first eight steps the bit length goes
1200 → 1202 → 1201 → 1203 → 1202 → 1203 → 1202 → 1204: a clear upward drift
where a typical orbit drifts down. It takes 15,782 steps, ratio 13.1517.

So the governing quantity is how many consecutive odd steps the bit pattern
sustains — all-ones maximizes it, alternating 01 collapses it in one move.
That is the right intuition for the width dependence too, but it does not
finish the job: bucketing w = 8 patterns by popcount gives only a +0.177
correlation with the ratio (means rise from 7.26 at three 1-bits to 7.97 at
six, but with wide overlap). The sign is right and the effect is real;
popcount alone is too crude to predict which widths deviate.

This changes the earlier conclusion. Whether a family can be *sustained*
slow — pushed far enough above baseline to bear on convergence at all — is
still open, and the mechanism is no better understood here than in base 10.

### Which construction is more anomalous: binary, decisively

The two families were never compared on the same footing — the decimal work
measured `steps / decimal-digit` at fixed decimal length, the binary work
measured `steps / log2(n)` at fixed bit length. Redone with `log2(n)` as the
denominator throughout and 40 seeds per family:

| construction | mean | random baseline | difference | z |
|---|---|---|---|---|
| decimal, 4-digit block, 400 digits | 7.1281 | 7.2003 | −0.072 | −0.85 |
| decimal, 4-digit block, 800 digits | 6.9963 | 7.1570 | −0.161 | −2.56 |
| decimal, 4-digit block, 1200 digits | 7.0289 | 7.1970 | −0.168 | −3.28 |
| **binary, w = 10, 1200 bits** | **8.0064** | 7.2473 | **+0.759** | **+3.54** |
| **binary, w = 4, 1200 bits** | **5.7978** | 7.3019 | **−1.504** | **−12.87** |

The binary effect is an order of magnitude larger in the quantity that
matters — the deviation itself. Decimal chains shift the ratio by 0.07–0.17;
binary widths shift it by 0.76 (w = 10) and 1.50 (w = 4). The z-scores are
closer than the deviations because the binary families also have wider
spread, but on effect size binary wins outright.

The decimal effect is also fragile in a way the binary one is not. Holding
*bit* length fixed instead of decimal length collapses it to z = −0.64 —
statistically nothing — while the same change leaves the binary results
untouched. So part of what the decimal comparison was detecting is an
artifact of what "same size" was taken to mean, not a property of the
numbers. The binary construction imposes structure on the 2-adic digits the
map actually operates on; the decimal one does not, and only leaks a weak
signal through whatever bit structure `n -> n*10^4 + s` incidentally leaves.

### Where the 1997 chain sits in all this: nowhere

Worth stating plainly, since the binary sweep invites the question. 1997 is
`0b11111001101`, **11 bits**, not 9. And the decimal chain is not a binary
repetition at all: `n -> n*10000 + 1997` with `10000 = 2^4 · 625` produces no
periodic bit pattern, as inspecting successive terms confirms. The chain
belongs to a different family than anything in the width table.

Taking 1997 as an 11-bit *binary* pattern instead and repeating it gives
ratio 7.5077, z = +0.72 — unremarkable. Width 11 as a whole is flat
(z = −0.33). Width 9, which would be the relevant one if 1997 were a 9-bit
value, is genuinely fast (z = −4.35), but 1997 is not a 9-bit value.

## Why this is not a lead on the conjecture

Three independent reasons, any one of which is sufficient.

**It measures the wrong thing.** The conjecture asserts that every positive
integer reaches 1. This finding says some numbers that reach 1 do so slightly
faster. A counterexample would be a number that *never* reaches 1. No amount
of precision about the speed of convergent trajectories constrains the
existence of non-convergent ones.

**The set has density zero.** Numbers built from a repeated 4-digit block
number about `10^4` out of `10^d` at d digits — a fraction of `10^-1196` at
1200 digits. Properties of a measure-zero family say nothing about the
integers as a whole.

**Average-case results are known not to suffice.** Terras (1976) established
density results for stopping times; Tao (2019) proved that almost all orbits
attain almost bounded values. The conjecture remains open regardless. The
gap between "almost all" and "all" is exactly where the difficulty lives, and
statistics about typical trajectories do not cross it.

## Errors made during this analysis

Recorded because the same mistake recurred in three different disguises, each
time by comparing numbers whose aggregation ranges or sample sizes differed.

**1. Comparing a mixed-range mean to a fixed-size mean.** The 1997 log's
overall mean, 23.1416, was compared against the theoretical 24.0118 and the
gap called significant. But the log's mean pools digit counts from 4 to
104,240, over half of its 26,060 entries above 50,000 digits, while the
theoretical value describes a fixed size. Restricting the log to comparable
ranges gives 23.59 (100–1,600 digits), 23.19 (1,600–10,000), 23.13
(50,000+) — the "constant" 23.13 is an artifact of where the samples are.

**2. Using the wrong sigma.** The gap was reported as 11σ using σ = 0.08,
the spread of the pooled 26,060-entry mean. The relevant quantity is the
spread of a *single* number's ratio at fixed size, which is σ ≈ 0.6–2.6 —
ten to thirty times larger. At 200 digits, fifteen random numbers ranged from
20.05 to 28.19. Corrected, the same gap is under 1σ.

**3. Declaring "no difference" from five points.** A z-test across five
digit sizes gave t = +1.15 and the effect was called absent. Repeating with
matched digit counts and larger samples gave t = −4.50, with the sign
consistent at all eight sizes tested. The five-point test was underpowered.

The intermediate claim that the ratio was systematically low *because* of
mod-4 structure, and the claim that 1997 was the highest among tested seeds,
are both withdrawn — the first is refuted above, the second rested on 16
cherry-picked iterations (multiples of 25) rather than the full range.

Separately, an earlier note in `collatz_1997_chain.py` predicted the ratio
should be about 8.0. That figure came from a model that ignored the forced
division after each `3n+1`; the correct prediction is 24.01, which the data
matches.

## What would actually be worth doing

Not continuing the searches. In order of value:

1. **Explain the width dependence in binary.** Widths 3, 5, 8, 10, 12 fall
   slow; width 4 falls very fast; 6, 16, 20, 24 show nothing. This is a
   concrete, self-contained puzzle with a clean algebraic setting: a number
   built from a repeated `w`-bit pattern is `p * (2^{wk} − 1)/(2^w − 1)`, so
   `2^w − 1` is the natural object to look at. The obvious first guess was
   tested and **fails**: whether 3 divides `2^w − 1` does not track the
   observed behavior (w = 8, 10, 12 are all divisible and slow, but so are
   w = 16, 20, 24, which show nothing; w = 3 and 5 are not divisible and are
   slow). Whatever governs this is subtler. Start here rather than with the
   base-10 version, which lacked this structure entirely.
2. **Read the literature.** Lagarias's annotated bibliography, Terras (1976),
   Tao (2019). The recurring lesson of this repository is that plausible
   ideas here are decades old, and reading is faster than rediscovering.

The scheduled GitHub Actions searches contribute to neither. The random
search has covered `10^-1019` of its space and last set a record on
2026-08-03; the cycle hunt covers `10^-136` per year at its best target.
