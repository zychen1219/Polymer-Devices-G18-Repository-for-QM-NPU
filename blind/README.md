# Blind test run

Date: 2026-11-15, 12:44 (local)
Data: `data/raw/blind_G18` (the released blind set, copied here unedited)
Command:

```bash
python run_blind.py --data data/raw/blind_G18 --out blind
```

Total wall-clock time: 13.1 s. No code was written for this run; `run_blind.py`
only calls the same `src/l1.py`, `src/l2.py` and `src/l3.py` that produce
`results/results.json`. That is the point of the exercise: the pipeline has to
run on a data set it has not seen.

## What it found

| phase | result |
| --- | --- |
| L1 | C = 3.4478, n_F = 4.102e11 cm^-2, mu_F = 29 301 cm^2/Vs |
| L2 | hard gap, Eg = 3.00 kT (77 meV), r_gap/r_cone = 0.126, gap fraction 1.000 over 200 resamples |
| L2 drift | +96.1 meV (index 0) to -135.2 meV (last index) |
| L3 | a = 0.0271 against tau = 0.030, delta95 = 2.57e-06, FWHM 1.74 pitch, PASS |
| L3 features | +0.0271 at (1,3), -0.0141 at (9,0) |

The blind wafer passes, the group's own wafer fails. The difference is where the
non-uniformity sits: the regular set carries a global gradient that reaches 3.2
per cent at a corner, the blind set carries a compact hot region at
`(x, y) = (-1.5, -2.0)` pitch with a peak of 2.7 per cent and a cold corner of
1.4 per cent. Two different wafers, two different answers from the same code.

## Checks performed during the run

- `data/raw/blind_G18` has the same eight files and the same schema as the
  regular set; only `c`, `EF_op_eV` and `rho_ohm_sq` differ;
- all 100 map values are finite and non-zero;
- `results.json` keeps the same six L3 keys, so the file is schema-identical to
  the regular run;
- the two L2 gates still fire on the blind data (zero crossing at index 31, cone
  wing exponent 1.979 against the gate of 2.0).

## Note on `pipeline_commit`

The repository did not exist when the blind run was first made, so the initial
log recorded `UNCOMMITTED-LOCAL-RUN`. The pipeline is deterministic (fixed seeds
`20261115` for the bootstrap and the noise floor), so the run was repeated after
the commits and the value was replaced with
`2008bef947f8bbf792a317c3c180ab9d2e95b80f`. Nothing else changed: the numbers
above are reproduced exactly.
