# Phase check — L2

Date: 2026-11-15
Scope: zero-dimensional inversion, group G18, 300 K
Command: `python reproduce.py`
Result: exit 0, all assertions passed

## Gates

| gate | what it proves | evidence |
| --- | --- | --- |
| neutral point before normalisation | the zero crossing is a property of the data, not of the scale | `R_H` changes sign between index 24 and 25 |
| sign change survives normalisation | the two curves were normalised separately | normalised `R_H` still has both signs |
| wing scaling of the cone | the model curve has the shape the guide requires | measured exponent 1.979 against the gate of 2.0 |
| reachability floor | a gap must beat the cone by 15 per cent | ratio 0.142 against 0.85 |
| gap floor | a gap must be at least one thermal energy | 4.25 kT against 1 kT |

The wing gate deserves a note. The first version of this check assumed
`R_H = (n - p) R_s^2` with `n - p` proportional to `n + p`, which predicts a
linear exponent. The measurement returned 1.98, the assumption was wrong, and
the gate now asserts the squared law that the guide actually states. The episode
is recorded in `docs/physics-notes.md` section 11 and in the AI error ledger.

## Result

```
model            hard gap
Eg               4.25 kT  (= 110 meV)
16-84 per cent   [4.25, 4.50] kT
r_cone           0.1116
r_gap            0.0159
r_gap / r_cone   0.142
gap fraction     1.000 over 200 resamples (seed 20261115)
E_F drift        +73.0 meV (index 0)  ->  -124.2 meV (index 55)
```

## Figures

`figures/l2_model_space.png`

- both model curves drawn in the normalised (R_xx, R_H) plane.
- the 56 group points lie on the red hard-gap curve and inside the grey cone.
- the visible separation is the evidence: the cone misses the low-resistance
  branch by about 0.05 in normalised units, the gap misses by about 0.01.

`figures/l2_trajectory.png`

- `E_F` against measurement index, bootstrap band shaded.
- the drift falls from +73 to -124 meV and crosses zero between indices 24 and 26;
  the small back-steps are the 1.1 meV quantisation of the E_F grid, not structure.
- the band is narrow at the start and widens to [-142, -120] meV at the end.

`figures/l2_repeats.png`

- the three raw repeats of both curves, with the sign-change window shaded.
- the repeats agree within their own scatter, so the neutral point is not an
  artefact of averaging.

## Gate sensitivity

`tests/test_l2_gates.py`, five cases, `python -m unittest discover -s tests`:

| mutation | expected | observed |
| --- | --- | --- |
| normalise both curves by one shared maximum | the residual must move | best-gap residual changes by more than 10 per cent |
| fold `R_H` to `|R_H|` | the neutral crossing must vanish | residual improves to 0.0148, trajectory ends at +128 meV, no zero crossing |
| data lying on the cone | the 15 per cent rule must keep the cone | model = cone, ratio >= 0.85 |
| `(n - p)` at two known points | sign follows the course convention | positive then negative |
| unmodified pipeline | all gates pass | 9 tests OK |

The second row is the interesting one: discarding the sign does not degrade the
fit at all, it slightly improves it, because the loop folds onto itself. What it
destroys is the assignment, and only the trajectory catches it. A pipeline that
checked residuals alone would have shipped that mistake.

## One thing the data shows that the pipeline does not need

`R_xx` peaks at index 20 while `R_H` crosses zero at index 24. The neutral point
is therefore defined by the Hall coefficient, not by the resistance maximum, and
the four-index separation is a property of this wafer rather than of the model.
The report should say so, because a reader who locates neutrality by the
resistance peak would place it four points early.

## Open items

1. The gap interval is one grid step wide. The width reflects the 0.25 kT grid,
   not the noise; the report must not present it as a confidence interval.
2. The model comparison is a two-hypothesis residual comparison. No likelihood
   ratio or information criterion is quoted, and none should be invented later.
3. `results/l2_diagnostics.json` carries the residuals and the bootstrap
   settings. It is not part of the schema in `results_template.json` and must
   never be merged into `results.json`.
4. The September checkpoint reported `1.25 kT`. The correction is described once
   in the L2 section of the report (see `docs/physics-notes.md` section 14).

## How to check this independently

Without running the code: take the two model curves at `Eg = 4.25 kT`, normalise
each by its own maximum, and compute the RMS distance from the 56 normalised
data points to the nearest model point. The gap curve must give `0.0159` and the
cone `0.1116`. Any residual below `0.01` or above `0.03` for the gap means the
grid, the normalisation or the sign convention has drifted from this record.
