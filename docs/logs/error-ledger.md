# AI error ledger — QXU6028 group G18

Each row records an error that came from an assistant draft, the check that exposed it, and the
correction. The last column gives the commit that carries each fix and the file it touches; `git show <hash> -- <file>`
shows the change, and `07a192d` carries the assertions that caught the first, third and fourth entries.

| # | Date / phase | Tool and task | AI claim or output | Evidence that exposed it | Correction and physical reason | Commit |
|---|---|---|---|---|---|---|
| 1 | 5 Nov, L1 | Assistant asked to derive the band-edge count | Reported C = 0.4 at E_F = 0.10 eV and concluded the brief's "C around 2" was wrong | Closed form evaluated at 73 meV returns 2.005; the full quadrature from 0 returns 2.599 at 100 meV | The two counts share one integrand and differ only in the lower limit (0 versus E_F). The draft had replaced the integral by D(E_F)k_BT, which truncates at the Fermi level in the wrong direction and gives a ratio that falls with bias instead of rising | `f2dae2a`, `src/l1.py`; gate in `07a192d` |
| 2 | 9 Nov, L1 | Assistant asked to make the third check work when no table is supplied | `ValueError: zero-size array to reduction operation maximum which has no identity` | Running `python reproduce.py` on a folder without `L1_data.csv` | The check now prints `[SKIP] gate 3: no measured table supplied` and the suite asserts that an empty table is a skip, not an exception. The physical content is unchanged: a check that cannot run must say so rather than crash | `f2dae2a`, `src/l1.py` (gate 3) |
| 3 | 12 Nov, L2 | Assistant asked to write the Hall coefficient of a two-carrier film | `R_H = (p − n)/(n + p)²` | Zero-crossing gate: the trajectory ran from −72 to +128 meV, so the first datum (R_H = +5.35e-12) was assigned a negative Fermi level on a wafer that starts n-type. The residual gate did not catch it: using the wrong sign lowered the residual to 0.0148 | Reversed to `R_H = (n − p)/(n + p)²`, which makes a positive Hall coefficient mean electron-like transport; the trajectory then runs from +73.0 to −124.2 meV and crosses zero at index 24 | `f2dae2a`, `src/l2.py`; gate in `07a192d` |
| 4 | 14 Nov, L3 | Assistant asked to apply the zero-mean constraint before solving | Removed the mean of `d`, then solved for the map | Residual rose from 1.8×10⁻⁴ to 1.1×10⁻³ (against a noise floor of 1.15×10⁻⁵) and the reported peak doubled, from 0.032 to 0.087 | The constraint acts on `m`. Row sums of A vanish, so a constant added to the map is invisible to the measurement; column sums do not vanish, so `A m` does reproduce the level of `d`. Centring the data deletes a signal the forward model can fit and forces the fit to spend its freedom on noise | `f2dae2a`, `src/l3.py`; test in `07a192d` |
| 5 | 14 Nov, L3 | Assistant asked to choose eta from the L-curve | Corner of maximum discrete curvature, at 62 eta0 | Residual at that eta is three orders of magnitude above the noise floor, and the roughness varies by only 15 per cent across four decades, so there is no interior corner to find | Reverted to the guide's reference eta0 = 1e-2‖A‖₂²/‖L‖₂², with the discrepancy-principle value 0.032 eta0 kept as a diagnostic. The reference value is what the guide specifies; the corner rule was an invention | `f2dae2a`, `src/l3.py` |

Two notes on what these entries show.

The first and third were not caught by the quantity the code was optimising. A wrong integration
limit produced a smooth, plausible curve, and a wrong sign of R_H produced a better residual.
Both were caught by a check on a property the fit does not see: the reference values of C, and the
sign of the first datum. That is why the zero-crossing gate is a gate and not a diagnostic.

The second and fourth are boundary and constraint errors. Neither changes the model; both change
what the numbers mean, and both were found by running the pipeline at its edges: with no table,
and with the constraint applied to the wrong variable.
