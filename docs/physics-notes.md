# Physics notes — L1 forward model

Working notes for the group. Every equation here is implemented in `src/l1.py`
and policed by an assertion; numbers quoted are from `results/results.json`.

## 1. The two counts

Both models describe the same band-edge density of states

```
g(E) = c E        (E >= 0, Dirac point at E = 0)
```

and the same occupancy, the Fermi-Dirac function

```
f(E) = 1 / (1 + exp[(E - E_F)/k_B T]).
```

They differ in one place only: the lower limit of the energy integral.

```
n_band(E_F) = N_A c ∫_0^{E_max}      E f(E) dE
n_F   (E_F) = N_A c ∫_{E_F}^{E_max}  E f(E) dE
```

The band-edge count integrates every state above the Dirac point. The
Fermi-count integrates only the states above the Fermi level, so it discards
the occupied states between 0 and E_F. Hence

```
n_band >= n_F,        C := n_band / n_F >= 1.
```

`C = 1` only at `E_F = 0`, where the two limits coincide. Any implementation
that returns `C < 1` has the limits the wrong way round; that is the first gate.

## 2. Closed form used as an independent check

With `E_max -> infinity`, the Fermi integral has no elementary antiderivative but
its low-temperature expansion is exact to the order we need:

```
∫_{E_F}^∞ E f(E) dE = k_B T [ E_F ln2 + (π^2/12) k_B T ] + O(k_B T/E_F)^2
```

so

```
n_F = N_A c k_B T ( E_F ln2 + π^2 k_B T / 12 ).
```

At 300 K, `k_B T = 25.852 meV`. The expansion is used only to police the
numerical quadrature: the two agree to `1.5e-10` relative at 50 meV and
`4.2e-9` at the working point. The expansion is never the value that enters
`results.json`.

## 3. The C factor and the mobility

Sheet resistance fixes the product of the two quantities one wants to separate:

```
rho = 1 / (e n mu)        =>        mu = 1 / (e n rho).
```

With `rho` held constant and the same measured conductance, the two models give

```
mu_F = 1/(e n_F rho),      mu_B = 1/(e n_band rho)
mu_F / mu_B = n_band / n_F = C.
```

So the corrected mobility is larger than the band-edge mobility by exactly the
same factor that the Fermi count is smaller. This is a statement about
degeneracy, not about a better measurement: at `E_F = 0.176 eV` the Fermi level
sits `6.8 k_B T` above the Dirac point, only the states immediately above `E_F`
conduct, and using the whole band over-counts the carriers by 4.48.

Units: `n` in cm^-2 and `rho` in ohm/sq give `mu` in cm^2/Vs with no extra
factor of `10^4`. The `10^4` appears only when `n` is quoted in m^-2 or `mu` in
m^2/Vs; mixing them is the most common error in this calculation.

## 4. Numbers at the working point

| quantity | value |
| --- | --- |
| c | 0.03875489491553498 |
| T | 300 K |
| E_F,op | 0.17610827554768008 eV |
| rho | 527.5077673801122 ohm/sq |
| n_band | 2.4558e12 cm^-2 |
| n_F | 5.4793e11 cm^-2 |
| C | 4.4819 |
| mu_B | 4818.0 cm^2/Vs |
| mu_F | 21594.0 cm^2/Vs |

Reference values of C, independent of `c` and of the group:

| E_F (eV) | C |
| --- | --- |
| 0 | 1.000 |
| 0.073 | 2.005 |
| 0.100 | 2.599 |
| 0.200 | 5.104 |

`C = 2` is reached at 73 meV, `C = 4` at 157 meV. The working point at 176 meV
therefore lies above both markers, which is why the headline number is 4.48 and
not the 2 quoted in the briefing's illustration.

## 5. Where the remaining uncertainty sits

Three sources, in decreasing size:

1. **The choice of count itself.** This is not an error bar, it is the 4.48
   factor: any comparison of literature mobilities must state which count was
   used. Nothing in the data can decide this; the model choice does.
2. **Upper integration limit.** Truncating at 2.5 eV instead of infinity drops a
   tail of relative size `exp[-(2.5 - E_F)/k_B T] ≈ 8e-40` at the working point.
   It is not a numerical concern at 300 K.
3. **Rounding of the quoted inputs.** `c`, `rho` and `E_F,op` are given to full
   double precision in `L1_params.json`; the table itself is reproduced to
   `7.6e-07`, so the forward model is not the limiting factor.

The quadrature uses 200 001 points on `[0, 2.5] eV`, giving a step of 12.5 µeV,
about 2000 points per `k_B T`. Halving the step changes `C` in the seventh
significant figure; the grid is therefore converged, not tuned.

## 6. Sensor statement

For a sensor whose response scales as `S ∝ 1/n`, replacing the band-edge count
with the Fermi count multiplies the predicted response by `C = 4.48` at this
bias. A design that quotes sensitivity from `n_band` will therefore
under-predict the measured signal by that same factor, and the discrepancy grows
with doping: at 200 meV it is 5.10, at 73 meV only 2.00.

## 7. Note on the September checkpoint

The checkpoint submitted in the live session reported `C ≈ 4.5` and
`mu ≈ 2.16e4`, consistent with the numbers above. Its code, however, carried two
defects worth recording for the AI error ledger:

- `N_A` was missing from the prefactor and the resulting scale was recovered by
  fitting `DOS_SCALE` against `L1_data.csv`. The fitted value equals `N_A` to
  five significant figures, which is why the checkpoint numbers were right and
  the code was still wrong: a free parameter was absorbing a constant that the
  methods guide already fixes.
- A second number, `C = 30.21`, appears in the same document from an
  implementation of `g(E) = c (E - E_F)`, a DOS that moves with the Fermi level
  and has no physical basis here. It is not used in this repository.

The present code removes both by construction: `N_A` is a constant in
`src/common.py` and the DOS is `c E` with the Dirac point fixed at zero, so the
CSV is reproduced without any fitted parameter.

## 8. L2 — the model pair and what can be identified

With mobility set to one and the carrier charge set to one, the two measured
curves for a given density of states are

```
R_H(E_F) = (n - p) / (n + p)^2        signed, never |R_H|
R_s(E_F) = 1 / (n + p)
n = N_A int_{E >= E_F} g(E) f(E) dE
p = N_A int_{E <= E_F} g(E) [1 - f(E)] dE
```

Any common factor in the mobility multiplies `n` and `p` alike and therefore
scales `R_s` and `R_H` by constants that the two separate normalisations remove.
This is why the guide allows `e = mu = 1`: the fit constrains the shape of the
pair, not its absolute size. The consequence is honest and worth stating once:
the inversion returns `E_F` and the gap, and it cannot return a mobility.

The two densities of states are

```
cone : g(E) = c |E|
gap  : g(E) = c (|E| - Eg/2)   for |E| > Eg/2,  0 otherwise
```

and each model curve is normalised by its own `max|R_H|` and its own `max R_s`.

## 9. The sign convention is fixed by the data, not by taste

The course convention is `(n - p)`: electrons are counted above `E_F` and holes
below it, so `R_H > 0` means `E_F > 0`. The first measured index has
`R_H = +5.35e-12` and the curve crosses zero between index 24 and 25, so the
sweep starts on the electron side and ends on the hole side.

Flipping to `(p - n)` mirrors the model curve about `R_H = 0` and, on this data,
returns a drift from `-72` to `+128 meV` — an electron-deficient start, which
contradicts the measured positive coefficient at index 0. The sign is therefore
not a free choice; it is pinned by the first data point and by the direction of
the drift.

## 10. Normalisation order

The neutral point is located first, from the sign change of the raw `R_H`
(index 24), and only then are the two curves normalised, separately. The order
matters: normalising before locating the zero lets one large outlier set the
scale and move the apparent crossing.

A test enforces the separation. Normalising both curves by a single shared
maximum changes the best-gap residual by more than 10 per cent, so the mistake
is detectable in the result, not only in the code.

## 11. Model comparison

Residuals on the group data, 56 indices averaged over three repeats:

| model | rms in normalised pair space |
| --- | --- |
| gapless cone | 0.1116 |
| hard gap, 1.25 kT | 0.0859 |
| hard gap, 4.25 kT | **0.0159** |

The minimum is interior, not at the edge of the scanned range: the residual
falls from 0.1034 at a half-gap of 0.25 kT, bottoms at 0.0159 (4.25 kT full
gap), and rises again to 0.0177 at 5 kT and 0.0208 at 8 kT. The ratio
`r_gap / r_cone = 0.142` is far below the 0.85 threshold, and `Eg = 4.25 kT`
is far above the 1 kT floor, so the hard gap is selected by both clauses of the
course rule.

Equivalently, `4.25 kT = 110 meV` at 300 K.

## 12. Bootstrap

200 resamples, seed 20261115, one of the three repeats drawn at random for each
index, the model library rebuilt once and reused. The hard gap is selected in
100 per cent of resamples, `Eg` has median 4.25 kT and a 16-84 per cent interval
of [4.25, 4.50] kT.

The interval is one grid step wide, and that is a real limitation rather than a
sign of precision: the gap grid is quantised in 0.25 kT steps, so the width of
the interval is set by the grid, not by the noise. Quoting more than two
significant figures for the gap would be false precision.

`E_F` drifts from `+73.0 meV` at index 0 to `-124.2 meV` at index 55, crossing
the Dirac point around index 25. The individual steps are quantised by the 400-point
E_F grid (1.1 meV), so repeated values and small back-steps are expected; the
trend, not the step sequence, carries the information. The bootstrap band is about `+/- 1
meV` at the start and `[-142, -120] meV` at the end: the fit is well constrained
near the gap edges and looser where the density of states is flat.

## 13. What this inversion cannot do

1. It compares two hypotheses. A smaller residual is evidence, not a
   likelihood ratio, and no goodness-of-fit statistic is quoted because the
   noise is not characterised well enough to justify one.
2. The 15 per cent rule is a course heuristic. A datum halfway between the two
   families would be reported as a cone, and the conservative direction is
   deliberate.
3. The residual curve is shallow near the minimum (0.0159 to 0.0177 between 4.25
   and 5 kT), so the gap is quoted with a width that reflects the grid rather
   than a confidence level.
4. Nothing here tests a gapped cone with a different exponent, or an asymmetric
   electron-hole density. Those models were not supplied by the course and are
   not compared.

## 14. Difference from the September checkpoint

The live-session submission reported a hard gap of `1.25 kT` with
`r_gap = 0.0207` and `82.5%` gap selection. Three things do not survive
re-computation:

- at `1.25 kT` the residual is `0.0859`, five times larger than at `4.25 kT`;
  the number `0.0207` does not correspond to that gap;
- the selection was not computed. The checkpoint code contains
  `flag = 'SELECTED (15% rule)' if (Eg == 1.25)` and the literal
  `0.020667 / r_cone` inside a print statement, so the label and the residual
  were typed in rather than derived;
- the reported residual ratio `0.185` is `0.0207 / 0.112`. Using the computed
  residual at that gap, the ratio is `0.77`, still under 0.85, so the September
  answer `gapped` is right even though its gap is wrong.

The present repository keeps the model class and the corrected gap, and reports
the correction once in the L2 section of the report rather than silently
changing the record.

## 15. L3 — what the forward model fixes

The supplied matrix `A` has 128 rows and 100 columns: one row per reading, one
column per map pixel. Each row has 50 non-zero entries, so a reading is a
weighted average over half the wafer rather than a point sample. The geometry is
fixed by `L3_meta.json`: a 10 x 10 grid over a region of 6 by 6 probe pitches,
64 placements of the probe, 2 orientations, and 128 = 64 x 2 readings.

Two properties of `A` decide the whole analysis. Its rows sum to zero to within
`1.7e-16`, so adding a constant to the map changes no prediction: the absolute
sheet conductivity is not in the data, only the shape of its spatial variation.
Its columns do not sum to zero, so the *level* of the readings is informative and
must be kept. These are different statements and conflating them costs a factor
of a hundred in fit quality.

## 16. The mean of `m`, not of `d`

The constraint `mean(m) = 0` follows from the row sums. It is enforced by
projecting the normal equations with `P = I - 11^T/N`, which is equivalent to
solving in the subspace orthogonal to the constant vector.

The trap is to read "zero mean" as a property of the data and to subtract
`mean(d)` before the fit. That is wrong. Because the columns of `A` do not sum
to zero, `A m` reproduces the level of `d`; removing that level leaves a target
that no admissible map can reach. With `d` centred the residual is `1.14e-3`,
against `1.15e-5` with the raw data, and the reconstructed peak inflates from
`0.032` to `0.087` as the fit chases noise. The check is now an assertion in
`tests/test_l3_gates.py`, so the regression cannot return quietly.

## 17. The regularisation parameter

The guide fixes a reproducible starting value,

    eta0 = 1e-2 ||A||_2^2 / ||L||_2^2 = 1.331e-3

with `L` the first-difference operator on the 10 x 10 grid and hence
`||L||_2 = 2.794`. It also asks for a sweep of `0.01, 0.1, 1, 10, 100` times
`eta0`, and for `eta0` to be used unless the sweep gives a documented reason to
change it.

Two candidate reasons were examined and neither survived.

The L-curve has no corner here. Over the full sweep the roughness `||Lm||` moves
only between `7.29e-2` and `6.16e-2`, a 15 per cent change, while the residual
moves across three decades. A discrete-curvature rule returns the end of the
scan, `62 eta0`, which is the flattest map rather than the corner.

The noise-level criterion is a real alternative: the repeated readings give
`sigma = 1.64e-6` per channel, so the mean of three carries `1.02e-6` and a
128-channel vector of that noise has norm `1.15e-5`. Matching it would need
`0.032 eta0`. But at `eta0` the residual is `1.79e-4`, only 15.5 times larger,
and the map differs by 1.2 per cent at the peak. The decision does not move.

Since the guide's default is inside the range where the answer is stable, it is
kept. The reported result is `eta = 1.00 eta0`, and the noise-level match is
recorded as a diagnostic in `results/l3_diagnostics.json` rather than used. This
is the conservative choice: a reader who prefers the noise-level value gets the
same map and the same verdict.

## 18. Resolution

The probe is a single-pixel pulse of unit amplitude, `q = delta - mean(delta)`,
pushed through the forward model and reconstructed with the same `eta` and the
same solver. The half-width of the recovered peak is the resolution, quoted in
probe pitches.

Three estimators are computed and only two are usable. Counting samples above
half maximum returns 1 pixel, which is an artefact: the peak falls between two
samples, so the count saturates at the grid. Linear interpolation of the two
crossings returns 1.67 pixels. A Gaussian fitted to the log of the three samples
around the peak returns 1.74. The quoted value is the larger of the two, 1.74,
which is the conservative reading.

The centre and off-centre probes agree to two digits, so the resolution is not a
strong function of position over the region. The ringing beside the peak, at
about -0.2 of its height, is the price of the first-difference penalty and sits
below the detection floor.

## 19. Detection floor

Every pixel has its own noise level, so the floor is found by simulation rather
than by a single number. Five hundred vectors of Gaussian noise with the measured
per-channel standard deviation are pushed through the same inverse operator, and
the 95th percentile of the per-map maximum is taken as `delta95 = 3.01e-05`
(seed 20261115). Because the statistic is a maximum over pixels, it already
contains the extra amplification that the corners receive; a per-pixel estimate
would understate them.

`delta95` is 0.09 per cent of the reference conductivity, two orders of magnitude
below the threshold, so the decision is not limited by noise.

## 20. The decision and its margin

The rule is the one in the guide: with `a = max|m|` over the zero-mean map, the
wafer passes if `a + delta95 < tau`, fails if `a - delta95 > tau`, and is
borderline otherwise, with `tau = 0.030`.

Here `a = 0.0324`, so `a - delta95 = 0.0323 > 0.030`: the wafer fails, and the
margin is 0.0024, about 8 per cent of the threshold. The honest description is
that the decision is statistically secure and practically narrow. The excess is
carried by the two corner pixels, at `+0.0324` and `-0.0311`, which are the least
constrained points of the map; the largest interior value is `0.0269`, below the
threshold. A reader who restricted the rule to the interior would pass the wafer.
The reconstruction says the non-uniformity is a single diagonal gradient of about
6.4 per cent peak to peak rather than a set of localised defects, and the report
says exactly that.
