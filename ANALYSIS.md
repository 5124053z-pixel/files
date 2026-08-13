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
| Is the effect real? | Yes, z = −3.52 at 1200 digits with matched sample sizes |
| Is it specific to 1997? | No. It appears in repeated-digit numbers generally |
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

No slow-falling family was found. Note that none of the tested block widths
came out *above* the random baseline either — every significant deviation in
the table above is negative. Whether a construction with r < 2 sustained over
a whole orbit exists is open, but it would have to constrain the trajectory
rather than just the starting value, and that is a much harder thing to
arrange.

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

1. **Explain the elevated r.** The mechanism behind repeated-digit numbers
   having r = 2.013 is open, and the non-monotonicity in block width (strong
   at 1, 3, 4; absent at 2, 5, 6) is a concrete, self-contained puzzle.
   Likely tractable with congruence arguments; would need the 2-adic
   structure of `x -> x*10^w + s` examined properly.
2. **Read the literature.** Lagarias's annotated bibliography, Terras (1976),
   Tao (2019). The recurring lesson of this repository is that plausible
   ideas here are decades old, and reading is faster than rediscovering.

The scheduled GitHub Actions searches contribute to neither. The random
search has covered `10^-1019` of its space and last set a record on
2026-08-03; the cycle hunt covers `10^-136` per year at its best target.
