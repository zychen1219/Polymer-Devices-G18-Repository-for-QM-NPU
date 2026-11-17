"""L2 zero-dimensional inversion: gapless cone versus symmetric hard gap.

Model (methods guide, section L2).  With e = mu = 1 the sheet resistance scales
out of the problem once R_H and R_s are normalised separately, so the two
measurable curves are

    R_H(E_F) = (n - p) / (n + p)^2          [signed, never abs]
    R_s(E_F) = 1 / (n + p)

with the Fermi-window populations

    n = N_A int_{E >= E_F} g(E) f(E) dE
    p = N_A int_{E <= E_F} g(E) [1 - f(E)] dE

and two candidate densities of states

    cone : g(E) = c |E|
    gap  : g(E) = c (|E| - Eg/2)  for |E| > Eg/2, else 0

Each model curve is normalised by its own max|R_H| and max R_s, the data pair by
the same rule.  Model selection uses the course rule: accept the hard gap only
when r_gap < 0.85 r_cone AND Eg >= 1 kT.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common import E_CHARGE, KB_EV, N_A, trapz

E_HALF = 2.5          # eV, symmetric DOS window
N_E = 4000
N_EF = 400
GAP_GRID_KT = np.round(np.arange(0.50, 8.00 + 1e-9, 0.25), 2)
BOOTSTRAP_N = 200
BOOTSTRAP_SEED = 20261115


# --------------------------------------------------------------------------
# raw data
# --------------------------------------------------------------------------
def load_repeats(path: str) -> dict:
    """Return (n_index, 3) arrays of R_H and R_xx plus the index column."""
    d = np.genfromtxt(path, delimiter=",", names=True)
    idx = d["measurement_index"].astype(int)
    n_idx = int(idx.max()) + 1
    reps = np.unique(d["repeat"])
    RH = np.full((n_idx, len(reps)), np.nan)
    RS = np.full((n_idx, len(reps)), np.nan)
    for k, rep in enumerate(reps):
        sel = d["repeat"] == rep
        RH[idx[sel], k] = d["R_H_raw"][sel]
        RS[idx[sel], k] = d["R_xx_raw"][sel]
    assert np.isfinite(RH).all(), "L2_data.csv has a missing or duplicated repeat"
    assert np.isfinite(RS).all(), "L2_data.csv has a missing or duplicated repeat"
    return dict(index=np.arange(n_idx), repeats=reps, RH=RH, RS=RS)


def normalise_pair(rh: np.ndarray, rs: np.ndarray) -> np.ndarray:
    """Normalise each curve by its own maximum absolute value (never together)."""
    rh = np.asarray(rh, float)
    rs = np.asarray(rs, float)
    out = np.column_stack([rh / np.max(np.abs(rh)), rs / np.max(rs)])
    assert np.isfinite(out).all(), "normalisation produced non-finite values"
    return out


# --------------------------------------------------------------------------
# model library
# --------------------------------------------------------------------------
def dos(E: np.ndarray, c: float, model: str, eg_kt: float, T: float) -> np.ndarray:
    if model == "cone":
        return c * np.abs(E)
    half = 0.5 * eg_kt * KB_EV * T
    return np.where(np.abs(E) > half, c * (np.abs(E) - half), 0.0)


def populations(E: np.ndarray, EFG: np.ndarray, g: np.ndarray, T: float):
    """Vectorised n(E_F) and p(E_F) over the whole E_F grid."""
    kT = KB_EV * T
    F = 1.0 / (1.0 + np.exp(np.clip((E[:, None] - EFG[None, :]) / kT, -700, 700)))
    n = N_A * trapz(np.where(E[:, None] >= EFG[None, :], g[:, None] * F, 0.0), E, axis=0)
    p = N_A * trapz(np.where(E[:, None] <= EFG[None, :], g[:, None] * (1.0 - F), 0.0), E, axis=0)
    return n, p


def model_pair(n: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Normalised (R_H, R_s) curve, keeping the sign of R_H."""
    s = n + p
    s_safe = np.where(s > 0, s, np.nan)
    rh = (n - p) / s_safe ** 2
    rs = 1.0 / s_safe
    keep = np.isfinite(rh) & np.isfinite(rs)
    rh, rs = rh[keep], rs[keep]
    return np.column_stack([rh / np.max(np.abs(rh)), rs / np.max(rs)])


def build_library(c: float, T: float) -> dict:
    E = np.linspace(-E_HALF, E_HALF, N_E)
    kT = KB_EV * T
    EFG = np.linspace(-8 * kT, 9 * kT, N_EF)
    gaps = np.asarray(GAP_GRID_KT, dtype=float)
    lib = {"EFG": EFG, "E": E, "T": T, "c": c, "gaps": gaps}
    for model, eg in [("cone", 0.0)] + [("gap", g) for g in gaps]:
        n, p = populations(E, EFG, dos(E, c, model, eg, T), T)
        lib[(model, eg)] = model_pair(n, p)
    return lib


# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------
def rms_to_curve(data_xy: np.ndarray, curve_xy: np.ndarray):
    d2 = ((data_xy[:, None, :] - curve_xy[None, :, :]) ** 2).sum(axis=2)
    return float(np.sqrt(np.mean(np.min(d2, axis=1)))), np.argmin(d2, axis=1)


def fit(data_xy: np.ndarray, lib: dict) -> dict:
    r_cone, j_cone = rms_to_curve(data_xy, lib[("cone", 0.0)])
    best = None
    for eg in lib["gaps"]:
        r, j = rms_to_curve(data_xy, lib[("gap", eg)])
        if best is None or r < best[1]:
            best = (eg, r, j)
    eg, r_gap, j_gap = best
    ratio = r_gap / r_cone
    use_gap = bool(ratio < 0.85 and eg >= 1.0)
    return dict(model="gap" if use_gap else "cone", Eg_kT=float(eg) if use_gap else 0.0,
                r_cone=r_cone, r_gap=r_gap, ratio=ratio,
                nearest=(j_gap if use_gap else j_cone),
                nearest_cone=j_cone, nearest_gap=j_gap)


def trajectory(nearest: np.ndarray, EFG: np.ndarray) -> np.ndarray:
    return EFG[nearest].copy()


# --------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------
def gate_zero_crossing(raw: dict) -> None:
    """Normalisation is per-curve, so the sign change must survive it."""
    rh = raw["RH"].mean(axis=1)
    s = np.sign(rh)
    flips = np.flatnonzero(s[:-1] * s[1:] < 0)
    print(f"  [PASS] R_H crosses zero between index {flips.tolist()} (sign change preserved by normalisation)")
    assert flips.size >= 1, "R_H never changes sign: the neutral point is undefined"
    xy = normalise_pair(rh, raw["RS"].mean(axis=1))
    assert np.any(xy[:, 0] > 0) and np.any(xy[:, 0] < 0), "normalised R_H lost its sign change"


def gate_wing_scaling(lib: dict, data_xy: np.ndarray, raw: dict) -> dict:
    """Gate from the methods guide: in the wings the ideal cone obeys
    |R_H~| proportional to R_xx~^2 (gamma = 2).

    The exponent is measured, not assumed.  A first attempt at this gate used the
    unipolar hand-wave R_H = (n - p) R_s^2 with n - p proportional to s, which
    predicts gamma = 1; the measurement rejected it (gamma = 1.98).  The reason is
    that over the populated part of the band n - p is far flatter than n + p, so
    the squared scaling of the guide is the correct one.
    """
    cone = lib[("cone", 0.0)]
    kT = KB_EV * lib["T"]
    efg = lib["EFG"]
    wing = np.abs(efg) > 4 * kT
    ok = np.isfinite(cone[:, 0]) & (np.abs(cone[:, 0]) > 1e-6) & (cone[:, 1] > 1e-6)
    sel = ok[wing]
    g_model = np.polyfit(np.log(cone[wing, 1][sel]), np.log(np.abs(cone[wing, 0][sel])), 1)[0]

    rs_d = data_xy[:, 1]
    rh_d = np.abs(data_xy[:, 0])
    keep = (rh_d > 0.2) & (rs_d > 0.2)
    g_data = np.polyfit(np.log(rs_d[keep]), np.log(rh_d[keep]), 1)[0]
    print(f"  [INFO] wing exponent gamma: model {g_model:.3f}, data {g_data:.3f} "
          f"(guide gates the cone at 2.0)")
    assert 1.5 < g_model < 2.5, "cone model does not follow |R_H| ~ R_xx^2 in the wings"
    return dict(gamma_model=float(g_model), gamma_data=float(g_data))


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def run(root: str, bootstrap_n: int = BOOTSTRAP_N, seed: int = BOOTSTRAP_SEED) -> dict:
    params = json.load(open(os.path.join(root, "L1_params.json"), encoding="utf-8"))
    c, T = float(params["c"]), float(params["T"])
    raw = load_repeats(os.path.join(root, "L2_data.csv"))

    print("L2 gate 1: neutral point survives normalisation")
    gate_zero_crossing(raw)

    lib = build_library(c, T)
    mu_rh = raw["RH"].mean(axis=1)
    mu_rs = raw["RS"].mean(axis=1)
    data_xy = normalise_pair(mu_rh, mu_rs)
    print("L2 gate 2: wing scaling of the calibrated cone")
    gamma = gate_wing_scaling(lib, data_xy, raw)

    centre = fit(data_xy, lib)
    ef_centre = trajectory(centre["nearest"], lib["EFG"])
    print(f"  r_cone = {centre['r_cone']:.5f}   r_gap = {centre['r_gap']:.5f} "
          f"({centre['Eg_kT']:.2f} kT)   ratio = {centre['ratio']:.3f}")
    print(f"  selected: {centre['model']}   Eg = {centre['Eg_kT']:.2f} kT")

    rng = np.random.default_rng(seed)
    n_idx, n_rep = raw["RH"].shape
    eg_samples, model_is_gap, ef_samples = [], [], []
    for _ in range(bootstrap_n):
        pick = rng.integers(0, n_rep, size=n_idx)
        rows = np.arange(n_idx)
        b_xy = normalise_pair(raw["RH"][rows, pick], raw["RS"][rows, pick])
        b = fit(b_xy, lib)
        eg_samples.append(b["Eg_kT"])
        model_is_gap.append(b["model"] == "gap")
        ef_samples.append(trajectory(b["nearest"], lib["EFG"]))
    eg = np.asarray(eg_samples)
    gap_frac = float(np.mean(model_is_gap))
    chosen = eg[eg > 0]
    if chosen.size:
        eg16, eg50, eg84 = np.percentile(chosen, [16, 50, 84])
    else:
        eg16 = eg50 = eg84 = 0.0
    ef_samples = np.asarray(ef_samples)
    ef16, ef50, ef84 = np.percentile(ef_samples, [16, 50, 84], axis=0)

    print(f"  bootstrap n={bootstrap_n}: gap fraction = {gap_frac:.3f}")
    print(f"  Eg median {eg50:.2f} kT, 16-84% [{eg16:.2f}, {eg84:.2f}] kT")
    print(f"  E_F trajectory {ef50[0]*1e3:+.1f} -> {ef50[-1]*1e3:+.1f} meV "
          f"(band at ends: [{ef16[0]*1e3:+.1f}, {ef84[0]*1e3:+.1f}] / "
          f"[{ef16[-1]*1e3:+.1f}, {ef84[-1]*1e3:+.1f}])")

    return dict(model=centre["model"], Eg_kT=float(eg50), Eg_uncertainty_kT=[float(eg16), float(eg84)],
                EF_trajectory_eV=[float(x) for x in ef_centre], gap_fraction=gap_frac,
                r_cone=centre["r_cone"], r_gap=centre["r_gap"], ratio=centre["ratio"],
                gamma=gamma, bootstrap_n=int(bootstrap_n), seed=int(seed),
                EF_band16_eV=[float(x) for x in ef16], EF_band84_eV=[float(x) for x in ef84],
                neutral_index=int(np.flatnonzero(np.sign(mu_rh[:-1]) * np.sign(mu_rh[1:]) < 0)[0]))


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def figures(root: str, res: dict, outdir: str) -> None:
    import os
    import matplotlib
    os.environ.setdefault("MPLCONFIGDIR",
                          os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".mplconfig"))
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    params = json.load(open(os.path.join(root, "L1_params.json"), encoding="utf-8"))
    c, T = float(params["c"]), float(params["T"])
    raw = load_repeats(os.path.join(root, "L2_data.csv"))
    lib = build_library(c, T)
    mu_rh, mu_rs = raw["RH"].mean(axis=1), raw["RS"].mean(axis=1)
    data_xy = normalise_pair(mu_rh, mu_rs)
    eg = res["Eg_kT"]
    sel = lib[("gap", eg)] if eg > 0 else lib[("cone", 0.0)]

    # 1. model space
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.plot(lib[("cone", 0.0)][:, 1], lib[("cone", 0.0)][:, 0], color="0.55", lw=1.2,
            label="gapless cone")
    if eg > 0:
        ax.plot(sel[:, 1], sel[:, 0], color="#c00000", lw=1.4, label=f"hard gap, {eg:.2f} kT")
    ax.plot(data_xy[:, 1], data_xy[:, 0], "o", ms=3.6, color="#1f4e79",
            label="G18 (mean of 3 repeats)")
    ax.axhline(0, color="0.8", lw=0.8)
    ax.set_xlabel(r"$\tilde R_{xx}$")
    ax.set_ylabel(r"$\tilde R_H$")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3, lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "l2_model_space.png"), dpi=300)
    plt.close(fig)

    # 2. Fermi trajectory
    ef = np.asarray(res["EF_trajectory_eV"]) * 1e3
    lo = np.asarray(res["EF_band16_eV"]) * 1e3
    hi = np.asarray(res["EF_band84_eV"]) * 1e3
    x = raw["index"]
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.fill_between(x, lo, hi, color="#1f4e79", alpha=0.18, lw=0, label="bootstrap 16-84%")
    ax.plot(x, ef, color="#1f4e79", lw=1.8, label="best fit")
    ax.axhline(0, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("measurement index")
    ax.set_ylabel(r"$E_F$ (meV)")
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.grid(alpha=0.3, lw=0.6)
    ax.set_title(f"{res['model']} fit, Eg = {eg:.2f} kT", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "l2_trajectory.png"), dpi=300)
    plt.close(fig)

    # 3. raw repeats around the sign change
    fig, axes = plt.subplots(2, 1, figsize=(6.2, 4.4), sharex=True)
    for k in range(raw["RH"].shape[1]):
        axes[0].plot(x, raw["RH"][:, k] * 1e12, lw=1.0, marker=".", ms=2.5, label=f"repeat {k + 1}")
        axes[1].plot(x, raw["RS"][:, k] * 1e12, lw=1.0, marker=".", ms=2.5)
    axes[0].axhline(0, color="0.4", lw=0.8)
    j = res["neutral_index"]
    axes[0].axvspan(j - 0.4, j + 0.4, color="#c00000", alpha=0.10, lw=0)
    axes[0].annotate("sign change", (j, 0), textcoords="offset points", xytext=(6, 8),
                     fontsize=7, color="#c00000")
    axes[0].set_ylabel(r"$R_H$ ($10^{-12}$ m$^2$/C)")
    axes[1].set_ylabel(r"$R_{xx}$ ($10^{-12}$ $\Omega$)")
    axes[1].set_xlabel("measurement index")
    axes[0].legend(fontsize=6.5, frameon=False, ncol=3, loc="upper right")
    for a in axes:
        a.grid(alpha=0.3, lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "l2_repeats.png"), dpi=300)
    plt.close(fig)
