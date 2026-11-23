"""Shared physics constants, loader and integration helpers for QXU6028 G18.

Every constant here is quoted by the course materials; nothing is fitted to the
data except where a comment says so explicitly.
"""
from __future__ import annotations

import json
import os

import numpy as np

# --- constants quoted by the methods guide -------------------------------
KB_EV = 8.617333262e-5          # eV/K
E_CHARGE = 1.602176634e-19      # C
N_A = 3.816e15                  # cm^-2, DOS prefactor fixed by the guide
E_MAX = 2.5                     # eV, upper integration limit
DEFAULT_T = 300.0

trapz = getattr(np, "trapezoid", np.trapz)


def fermi_dirac(E: np.ndarray, EF: float, T: float = DEFAULT_T) -> np.ndarray:
    """f(E) = 1/(1+exp((E-EF)/kT)); clipped so exp never overflows."""
    x = np.clip((E - EF) / (KB_EV * T), -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(x))


def grid(n: int = 200001) -> np.ndarray:
    return np.linspace(0.0, E_MAX, n)


def load_group(root: str) -> dict:
    """Load the seven per-group files from a directory."""
    def js(name):
        with open(os.path.join(root, name), encoding="utf-8") as fh:
            return json.load(fh)

    params = js("L1_params.json")
    l2meta = js("L2_meta.json")
    l3meta = js("L3_meta.json")

    l1 = np.genfromtxt(os.path.join(root, "L1_data.csv"), delimiter=",", names=True)
    l2 = np.genfromtxt(os.path.join(root, "L2_data.csv"), delimiter=",", names=True)
    A = np.load(os.path.join(root, "L3_A.npy"))
    d = np.loadtxt(os.path.join(root, "L3_d.csv"), delimiter=",")

    return dict(params=params, l2meta=l2meta, l3meta=l3meta,
                l1=l1, l2=l2, A=A, d=d)


def assert_close(value, target, tol, label):
    ok = abs(value - target) <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {value:.6g} vs {target:.6g} (tol {tol:g})")
    assert ok, f"{label} failed: {value} vs {target}"
    return ok


def style() -> None:
    """Print-legible defaults: the report embeds each figure near its native size."""
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 12,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "axes.linewidth": 0.8,
        "grid.alpha": 0.35,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
