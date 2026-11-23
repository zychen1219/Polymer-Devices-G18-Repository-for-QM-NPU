# QXU6028 G18 — carrier counting, gap inversion and wafer mapping

Coursework repository for QXU6028 (Polymer Devices), group G18:
Qiao Qinpeng, Wang Ziyue, Jiang Tianyu, Chen Zhengyi.

The repository holds one pipeline that is run on the group's synthetic
measurement set. Each phase replaces the physics in the middle of the pipeline
and keeps the same file layout, so `results.json` never changes shape.

```
raw data -> forward model -> inversion -> uncertainty -> decision -> report
              (L1)          (L2, L3)     (resampling)   (pass/fail)
```

## Run it

```bash
python -m pip install -r requirements.txt
python reproduce.py                 # writes results/ and figures/
python -m unittest discover -s tests -v   # proves the gates can fail
```

`reproduce.py` takes no arguments. It writes `results/results.json`, seven
figures and two diagnostics files, and it stops with a non-zero exit
code if any phase gate fails. Every printed `[PASS]` line is an assertion in `src/`, not a
comment.

Optional flags: `--data <dir>` to point at another group directory, `--out
<dir>` to write elsewhere.

## Layout

```
reproduce.py            one-command driver
src/common.py           constants (N_A, k_B, E_MAX), loader, integration helpers
src/l1.py               L1 forward model: two carrier counts and the C factor
src/l2.py               L2 inversion: cone versus hard gap, bootstrap, trajectory
src/l3.py               L3 inversion: regularised map, resolution, detection floor
data/raw/regular_G18/   group data exactly as distributed (never edited)
results/results.json    machine-readable results, written only by reproduce.py
results/l2_diagnostics.json   residuals and bootstrap settings, not part of the schema
results/l3_diagnostics.json   eta sweep, resolution probes, noise floor, not in the schema
figures/                figures, written only by reproduce.py
docs/physics-notes.md   model derivations and the reasoning behind each gate
reports/check-*.md      phase acceptance notes: what passed and what is open
docs/logs/              AI usage log, error ledger, self-assessment, contribution statement
tests/test_gates.py     L1: breaks the physics on purpose and requires the gates to fire
tests/test_l2_gates.py  L2: normalisation order, sign convention and the gap rule
tests/test_l3_gates.py  L3: row sums, the mean of d, eta stability, FWHM, noise floor
run_blind.py            same pipeline on the released blind set, with timings
tools/build_report_docx.py   rebuilds the report document from its markdown source
data/raw/blind_G18/     blind set exactly as released
blind/                  blind results, run log and its figures
```

The written report is submitted through QMplus and is deliberately not stored here. `tools/build_report_docx.py` is kept so that the document can be rebuilt from its markdown source, and every number quoted in it can be checked against `results/results.json`.

## Status

| phase | state | headline number |
| --- | --- | --- |
| L1 forward model | done | C = 4.4819 at E_F = 0.1761 eV |
| L2 zero-dimensional inversion | done | hard gap, Eg = 4.25 kT (110 meV) |
| L3 two-dimensional inversion | done | FAIL, a = 0.0324 against 0.030 |
| blind test | done | PASS, a = 0.0271 against 0.030 |
| report | done, submitted separately | 2287 words, 9 pages; the source and the PDF are not in this repository |

## Results (L1)

| quantity | value |
| --- | --- |
| c (DOS coefficient) | 0.0387549 |
| working point E_F | 0.176108 eV |
| n_band | 2.4558e12 cm^-2 |
| n_F | 5.4793e11 cm^-2 |
| C = n_band / n_F | 4.4819 |
| mu from Fermi count | 21594 cm^2/Vs |
| mu from band-edge count | 4818 cm^2/Vs |
| sheet resistance (fixed) | 527.51 ohm/sq |

The table in `data/raw/regular_G18/L1_data.csv` is reproduced from `c` and the
guide value of `N_A` alone, to a maximum relative error of 7.6e-07. Nothing is
fitted.

## Results (L2)

| quantity | value |
| --- | --- |
| model | hard gap |
| Eg | 4.25 kT = 110 meV |
| Eg bootstrap 16-84% | 4.25 - 4.50 kT |
| r_cone | 0.1116 |
| r_gap | 0.0159 |
| r_gap / r_cone | 0.142 |
| hard-gap selection rate | 200/200 resamples |
| E_F drift | +73.0 meV -> -124.2 meV |

The gap is selected by both clauses of the course rule: the residual ratio is
below 0.85 and the gap is above 1 kT. The 16-84 per cent interval is one grid
step wide because the gap grid is quantised in 0.25 kT steps; it is a resolution
limit, not a confidence statement.

The September checkpoint reported 1.25 kT for the same data. The residual at
that gap is 0.0859, five times the minimum, and the checkpoint's selection label
was a literal in its own source. The correction is described once in
`docs/physics-notes.md` section 14.

## Results (L3)

| quantity | value |
| --- | --- |
| eta | 1.00 eta0 = 1.331e-3 |
| residual | 1.79e-04 (15.5 x the noise level of the mean) |
| map range | -0.0311 at (0,0) to +0.0324 at (9,9) |
| peak-to-peak non-uniformity | 6.4 per cent |
| FWHM | 1.74 probe pitch |
| delta95 | 3.01e-05 (500 noise reconstructions, seed 20261115) |
| a = max abs map value | 0.0324 |
| decision | FAIL, confident, margin 0.0024 |

The reconstruction is a single diagonal gradient rather than a set of localised
defects. The verdict rests on the two corner pixels; the largest interior value
is 0.0269, below the threshold, and `reports/check-l3.md` says so explicitly.
`eta0` is the guide's reference value and is kept because the map moves by 1.2
per cent if the noise-level match is used instead, so the parameter is not
carrying the answer.

## Blind test

`run_blind.py` runs the same three modules on the blind set and writes
`blind/blind_results.json`, `blind/run_log.txt` and `blind/pipeline_commit.txt`.
Thirteen seconds end to end, no code changes. The blind wafer passes
(`a = 0.0271` against `0.030`) with a compact hot region at (-1.5, -2.0) pitch,
where the group's own wafer fails on a global gradient. `blind/README.md` has the
full comparison.
