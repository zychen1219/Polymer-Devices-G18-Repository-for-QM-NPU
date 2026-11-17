# Phase check — L3

Date: 2026-11-15
Scope: two-dimensional inversion, 10 x 10 map over the 6s x 6s region, group G18
Command: `python reproduce.py`
Result: exit 0, all assertions passed

## Gates

| gate | what it proves | evidence |
| --- | --- | --- |
| row sums of A | the forward model is blind to a constant in the map, so the mean is not identifiable | max abs row sum 1.67e-16 |
| shapes and finiteness | the data are the 3 x 128 readings the guide describes | `(3, 128)`, all finite, range -0.0119 to 0.0203 |
| zero mean of the solution | the constraint is enforced by the projector, not by a post-hoc shift | mean(m) = 0.0e+00 for two independent solves |
| mean of d is *not* removed | the level of the data belongs in the fit | residual 1.15e-05 with d, 1.14e-03 with d centred |
| stability of the verdict | the answer must not hinge on a tuned parameter | FAIL at 0.1, 1 and 10 eta0 |
| resolution beats the grid | the map is not a resampling of the 10 x 10 grid | FWHM 1.74 pitch, above 1 pixel |
| signal above the floor | the judged extreme is not a noise excursion | 0.0324 against delta95 = 3.0e-05 |

## Result

```
eta                1.00 eta0  (= 1.3310e-03)
residual          1.79e-04   (15.5 x the noise level of the averaged data)
noise level       1.15e-05   (sigma of the repeat mean, 128 channels)
map range         -0.0311 (pixel 0,0)  ->  +0.0324 (pixel 9,9)
FWHM              1.74 probe pitch (centre and off-centre probes agree)
delta95           3.01e-05 over 500 noise reconstructions (seed 20261115)
a = max|m|        0.0324
a - delta95       0.0323   >  tau = 0.030
decision          FAIL (confident)
```

`results.json` carries the 100 map values, the four features above the floor, the
resolution in probe pitches, the detection limit, and the boolean decision.

## Figures

`figures/l3_map.png`

- the reconstruction is a single diagonal ramp, not a set of isolated spots;
- it rises from -3.1 per cent at the bottom-left corner to +3.2 per cent at the
  top-right, so the peak-to-peak non-uniformity is about 6.4 per cent;
- the four features above the floor are the two corners (+0.0324 at 9,9 and
  -0.0311 at 0,0) and two adjacent points (+0.0269 at 7,9 and -0.0256 at 2,0),
  which are the shoulders of the same ramp rather than separate defects.

`figures/l3_lcurve.png`

- residual against roughness over 81 values of eta from 0.01 to 100 eta0;
- the dashed line is the noise level of the averaged data;
- the chosen point sits at 1.0 eta0, 15.5 times above that line;
- the curvature of this curve has no interior maximum: the roughness spans only
  6.2 to 7.3 (in units of 1e-2) while the residual spans three decades, so the
  discrete corner rule would return the end of the scan (62 eta0). This is why
  the guide's default eta0 is used and the noise-level match is reported as a
  diagnostic instead: it would give 0.032 eta0, move the map by 1.2 per cent at
  the peak, and leave the decision unchanged.

`figures/l3_resolution.png`

- the reconstructed response to a single-pixel pulse injected at the centre and
  off-centre, each normalised to its own peak;
- both profiles fall to half maximum between the sample either side of the peak,
  which is why counting samples gives 1 pixel and interpolation gives 1.67;
- the negative shoulders at -0.2 of the peak are the regularisation's ringing,
  not a second feature: they sit below the detection floor.

## Gate sensitivity

| eta / eta0 | residual | \|\|Lm\|\| | max\|m\| | decision |
| --- | --- | --- | --- | --- |
| 0.01 | 9.37e-06 | 7.29e-02 | 0.0328 | FAIL |
| 0.10 | 2.23e-05 | 7.28e-02 | 0.0327 | FAIL |
| 1.00 | 1.79e-04 | 7.25e-02 | 0.0324 | FAIL |
| 10.0 | 1.41e-03 | 7.05e-02 | 0.0312 | FAIL |
| 100 | 8.66e-03 | 6.16e-02 | 0.0271 | PASS |

The verdict flips only at 100 eta0, where the residual is 800 times the noise
level and the map has visibly flattened. That point is over-smoothed by the
standard the guide sets, so it is not used. The honest statement in the report
is narrower than "the wafer fails": the excess is 10 per cent above the threshold
and it is carried entirely by two corner pixels, while the interior maximum is
0.0269, below the threshold.

## Open items

- the corners are the least constrained pixels in the map; the noise propagation
  accounts for their larger gain through delta95, but the report should say so;
- the decision rule is applied to the peak of the zero-mean deviation map, which
  is what the guide specifies; a reader who applied it to the interior only would
  obtain PASS, and the report states this explicitly;
- no L-curve corner is claimed. The parameter is the guide's reference value and
  the sensitivity table is the evidence that it is not carrying the answer.

## How to check this independently

1. run `python reproduce.py` and read the L3 block it prints;
2. `python -m unittest discover -s tests -v` (18 tests, all pass);
3. change `TAU` in `src/l3.py` to 0.04 and rerun: the decision must become PASS
   while the map is unchanged, which shows the threshold is the only thing
   separating the two answers;
4. multiply the data by 1.01 and rerun: the map must scale by 1.01, which shows
   the operator is linear in d and the zero-mean step is not absorbing signal.
