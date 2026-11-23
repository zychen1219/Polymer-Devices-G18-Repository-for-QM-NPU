"""L1 forward model: band-edge count versus Fermi-referenced count.

Both counts share the same DOS g(E) = c*E and the same occupancy f(E).  They
differ only in the lower limit of integration:

    n_band(E_F) = N_A * c * int_0^{Emax}   E f(E) dE
    n_F   (E_F) = N_A * c * int_{E_F}^{Emax} E f(E) dE

so C = n_band / n_F is >= 1 and the Fermi-referenced count is the smaller one.
"""
from __future__ import annotations

import numpy as np

from common import style, E_CHARGE, E_MAX, KB_EV, N_A, DEFAULT_T, fermi_dirac, grid, trapz

_E = grid()
_CACHE: dict[tuple, float] = {}


def n_band(EF: float, c: float, T: float = DEFAULT_T) -> float:
    """Band-edge count: integrate from the Dirac point (E = 0)."""
    f = fermi_dirac(_E, EF, T)
    return N_A * c * trapz(_E * f, _E)


def n_fermi(EF: float, c: float, T: float = DEFAULT_T) -> float:
    """Fermi-referenced count: integrate from E_F upwards."""
    E = _E[_E >= EF]
    f = fermi_dirac(E, EF, T)
    return N_A * c * trapz(E * f, E)


def n_fermi_closed(EF: float, c: float, T: float = DEFAULT_T) -> float:
    """Independent analytic value of n_F, used only to police the integral.

    With the upper limit taken as infinity the integral is
        kT [ E_F ln2 + pi^2 kT / 12 ].
    """
    kT = KB_EV * T
    return N_A * c * kT * (EF * np.log(2.0) + np.pi ** 2 * kT / 12.0)


def correction(EF: float, c: float, T: float = DEFAULT_T) -> float:
    nb, nf = n_band(EF, c, T), n_fermi(EF, c, T)
    return nb / nf if nf > 0 else np.nan


def mobility(n_cm2: float, rho_ohm_sq: float) -> float:
    """mu = 1/(e n rho).  n in cm^-2 and rho in ohm/sq give cm^2/Vs directly."""
    return 1.0 / (E_CHARGE * n_cm2 * rho_ohm_sq)


def sweep(c: float, T: float, efs: np.ndarray) -> np.ndarray:
    return np.array([correction(ef, c, T) for ef in efs])


def run(params: dict, l1_table, verbose: bool = True) -> dict:
    c, T, rho = params["c"], float(params["T"]), float(params["rho_ohm_sq"])
    EF_op = float(params["EF_op_eV"])

    if verbose:
        print("L1 gate 1: reference values of C (independent of c)")
    for ef, target, tol in [(0.0, 1.000, 0.002), (0.073, 2.005, 0.010),
                            (0.10, 2.599, 0.010), (0.20, 5.104, 0.020)]:
        got = correction(ef, c, T)
        ok = abs(got - target) < tol
        print(f"  [{'PASS' if ok else 'FAIL'}] C({ef:.3f} eV) = {got:.4f} (guide {target:.3f})")
        assert ok, "reference C value failed"

    if verbose:
        print("L1 gate 2: numeric integral matches the closed form")
    for ef in (0.05, 0.176):
        a, b = n_fermi(ef, c, T), n_fermi_closed(ef, c, T)
        rel = abs(a - b) / b
        print(f"  [{'PASS' if rel < 5e-3 else 'FAIL'}] n_F({ef}) = {a:.4e}, closed {b:.4e}, rel {rel:.2e}")
        assert rel < 5e-3

    n_rows = 0 if l1_table is None else len(l1_table)
    if n_rows == 0:
        print("  [SKIP] gate 3: no measured table supplied")
    else:
        if verbose:
            print("L1 gate 3: table reproduced from c and N_A alone (no calibration)")
        rel = []
        for row in l1_table:
            rel.append(abs(n_band(float(row["EF_eV"]), c, T) - float(row["n_band_cm2"]))
                       / float(row["n_band_cm2"]))
        rel = np.array(rel)
        print(f"  [{'PASS' if rel.max() < 1e-3 else 'FAIL'}] max rel err over {n_rows} rows = {rel.max():.3e}")
        assert rel.max() < 1e-3

    C = correction(EF_op, c, T)
    nb = n_band(EF_op, c, T)
    nf = n_fermi(EF_op, c, T)
    mu_b = mobility(nb, rho)
    mu_f = mobility(nf, rho)
    if verbose:
        print(f"L1 operating point E_F = {EF_op:.4f} eV")
        print(f"  n_band = {nb:.4e} cm^-2   n_F = {nf:.4e} cm^-2")
        print(f"  C = {C:.4f}")
        print(f"  mu_band = {mu_b:.1f}   mu_F = {mu_f:.1f} cm^2/Vs   ratio = {mu_f / mu_b:.4f}")
    assert mu_f / mu_b == C or abs(mu_f / mu_b - C) < 1e-6, "mu ratio must equal C"

    return dict(correction_factor=float(C), n_fermi_op_cm2=float(nf),
                mobility_corrected_cm2_Vs=float(mu_f),
                n_band_op_cm2=float(nb), mobility_band_cm2_Vs=float(mu_b),
                rho_ohm_sq=rho, EF_op_eV=EF_op)


def figure(params: dict, out_png: str) -> None:
    import os
    import matplotlib
    # keep the font cache inside the project so nothing is written to the user profile
    os.environ.setdefault("MPLCONFIGDIR",
                          os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".mplconfig"))
    matplotlib.use("Agg")
    style()
    import matplotlib.pyplot as plt

    c, T = params["c"], float(params["T"])
    EF_op = float(params["EF_op_eV"])
    efs = np.linspace(0.0, 0.30, 301)
    C = sweep(c, T, efs)

    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(efs, C, color="#1f4e79", lw=2)
    for target, ls in ((2.0, ":"), (4.0, ":")):
        idx = int(np.argmin(np.abs(C - target)))
        ax.axhline(target, color="0.55", ls=ls, lw=0.9)
        ax.plot(efs[idx], C[idx], "o", ms=5, mfc="white", mec="#1f4e79")
        ax.annotate(f"C = {target:.0f} at {efs[idx] * 1e3:.0f} meV",
                    (efs[idx], C[idx]), textcoords="offset points", xytext=(8, -12),
                    fontsize=8, color="0.25")
    ax.plot(EF_op, correction(EF_op, c, T), "s", ms=6, color="#c00000")
    ax.annotate(f"working point\n{correction(EF_op, c, T):.2f}",
                (EF_op, correction(EF_op, c, T)), textcoords="offset points",
                xytext=(-10, 14), ha="right", fontsize=8, color="#c00000")
    ax.set_xlabel(r"$E_F$ (eV)")
    ax.set_ylabel(r"$C = n_{band}/n_{F}$")
    ax.set_xlim(0, 0.30)
    ax.grid(alpha=0.3, lw=0.6)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
