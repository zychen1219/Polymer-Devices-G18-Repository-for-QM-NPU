# Phase check — L1

Date: 2026-11-13
Scope: forward model, group G18, 300 K
Command: `python reproduce.py`
Result: exit 0, all assertions passed

## Gates

| gate | what it proves | evidence |
| --- | --- | --- |
| C(0) = 1.000 | the two integrals coincide at the Dirac point | 1.0000 vs guide 1.000 |
| C(0.073) = 2.005 | the C = 2 marker sits where the guide says | 2.0049 vs 2.005 |
| C(0.100) = 2.599 | the briefing's 0.10 eV example | 2.5988 vs 2.599 |
| C(0.200) = 5.104 | the curve keeps rising with doping | 5.1044 vs 5.104 |
| numeric vs closed form | the quadrature is right, not just self-consistent | rel 1.5e-10 (50 meV), 4.2e-9 (176 meV) |
| CSV reproduced | no free parameter is hiding in the model | max rel err 7.6e-07 over 12 rows |
| mu_F / mu_B = C | the mobility statement is arithmetically forced by fixed rho | 4.4819 / 4.4819 |

The first gate is the one the briefing warns about: if any of these came out
below 1, the integration limits would be transposed.

## Operating point

```
E_F   = 0.17610825 eV          (L1_params.json, not the nearest table row)
n_band= 2.4558e12 cm^-2
n_F   = 5.4793e11 cm^-2
C     = 4.4819
mu_F  = 21594.0 cm^2/Vs        mu_B = 4818.0 cm^2/Vs
```

## Figure

`figures/l1_correction.png`

- x axis 0 to 300 meV, y axis C from 1 upward.
- open circles mark C = 2 (73 meV) and C = 4 (157 meV).
- red square is the working point at 176 meV, C = 4.48.
- the curve is monotonic and starts at exactly 1, so the limit check is visible,
  not only asserted.

## Gate sensitivity

`tests/test_gates.py` checks that the gates can fail, not only that they pass.
Four cases, `python -m unittest discover -s tests -v`, all `ok`:

| mutation | expected | observed |
| --- | --- | --- |
| swap the two integration limits | C < 1, gate 1 raises | C(0.073) = 0.4988, `reference C value failed` |
| drop `N_A` from the prefactor | gate 3 raises | max rel err = 1.000, raises |
| unmodified model | all gates pass | 7 PASS lines, exit 0 |
| working point | `mu_F / mu_B = C` | ratio 4.4819 vs C 4.4819 |

The first case carries a lesson worth keeping: at `E_F = 0` the swapped limits
still give `C = 1.000`, because both integrals start and end at the same place.
A single reference point at the Dirac point would have passed the broken code.
The 73 meV point is what actually catches it, which is why the gate uses four
points and not one.

## Open items

1. The figure is a single panel. If the report needs a second panel, the natural
   one is `mu_F` versus `E_F` at fixed rho, which is just `C(E_F)` rescaled.
2. `pipeline_commit` reads `UNCOMMITTED-LOCAL-RUN` until the first commit is
   made. After the first push, re-run `reproduce.py` so the hash in
   `results.json` names the commit that produced the numbers.
3. The September checkpoint quoted `C ≈ 4.5`; the present value is 4.4819, so
   the two agree to two significant figures. The difference must still be
   described once in the report (see `docs/physics-notes.md` section 7).

## How to check this independently

Without running any code: take `n_F = N_A c k_B T (E_F ln2 + π^2 k_B T / 12)`,
substitute `N_A = 3.816e15`, `c = 0.0387549`, `k_B T = 0.025852 eV`,
`E_F = 0.176108`, and the result must be `5.479e11`. Then integrate `c E f(E)`
from 0 to 2.5 eV for `n_band` and divide. Any discrepancy above 1e-3 means the
code and this note disagree.
