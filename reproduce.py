#!/usr/bin/env python3
"""QXU6028 G18 - one-command reproduction.

    python reproduce.py

Writes results/results.json plus the figures from data/raw/regular_G18.
Exit code 0 means every phase gate passed; any failed assertion stops the run.

Optional:
    python reproduce.py --data <dir>   use a different data directory
    python reproduce.py --out  <dir>   write outputs somewhere else
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
# matplotlib and its wheels are vendored so a marker needs no internet
sys.path.insert(0, os.path.join(HERE, ".deps"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mplconfig"))

import l1  # noqa: E402
import l2  # noqa: E402
import l3  # noqa: E402
from common import load_group  # noqa: E402

DEFAULT_DATA = os.path.join(HERE, "data", "raw", "regular_G18")


def pipeline_commit() -> str:
    """Frozen hash of the commit that produced this run (brief section 7.3).

    Only accepted when this directory *is* the repository root. Otherwise the
    enclosing repository (an editor workspace, a copy) would be credited with a
    hash that does not exist in the coursework repository.
    """
    def git(*args):
        return subprocess.run(["git", *args], cwd=HERE,
                              capture_output=True, text=True, timeout=10)

    try:
        top = git("rev-parse", "--show-toplevel")
        if top.returncode != 0:
            return "UNCOMMITTED-LOCAL-RUN"
        if os.path.normcase(os.path.realpath(top.stdout.strip())) != \
                os.path.normcase(os.path.realpath(HERE)):
            return "UNCOMMITTED-LOCAL-RUN"
        head = git("rev-parse", "HEAD")
        if head.returncode == 0 and head.stdout.strip():
            return head.stdout.strip()
    except Exception:
        pass
    return "UNCOMMITTED-LOCAL-RUN"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--out", default=HERE)
    args = ap.parse_args()

    data = load_group(args.data)
    params = data["params"]
    print(f"QXU6028 {params['group']}  T = {params['T']} K  "
          f"E_F,op = {params['EF_op_eV']:.4f} eV")

    res = l1.run(params, data["l1"])
    res2 = l2.run(args.data)

    os.makedirs(os.path.join(args.out, "figures"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "results"), exist_ok=True)
    res3 = l3.run(args.data)
    l1.figure(params, os.path.join(args.out, "figures", "l1_correction.png"))
    l2.figures(args.data, res2, os.path.join(args.out, "figures"))
    l3.figures(args.data, res3, os.path.join(args.out, "figures"))

    with open(os.path.join(args.data, "results_template.json"), encoding="utf-8") as fh:
        results = json.load(fh)

    results["pipeline_commit"] = pipeline_commit()
    results["L1"] = {
        "correction_factor": round(res["correction_factor"], 4),
        "n_fermi_op_cm2": res["n_fermi_op_cm2"],
        "mobility_corrected_cm2_Vs": round(res["mobility_corrected_cm2_Vs"], 1),
    }
    eg16, eg84 = res2["Eg_uncertainty_kT"]
    results["L2"] = {
        "model": res2["model"],
        "Eg_kT": round(float(res2["Eg_kT"]), 2),
        "Eg_uncertainty_kT": round(float(eg84 - eg16) / 2.0, 2),
        "EF_trajectory_eV": [round(float(x), 6) for x in res2["EF_trajectory_eV"]],
    }
    results["L3"] = {
        "map": [round(float(v), 6) for v in res3["map"]],
        "features": [{k: (round(float(v), 6) if isinstance(v, float) else v)
                      for k, v in f.items()} for f in res3["features"]],
        "resolution_probe_pitch": round(float(res3["resolution_probe_pitch"]), 2),
        "detection_limit_fraction": round(float(res3["detection_limit_fraction"]), 8),
        "pass": res3["pass_"],
        "confidence": res3["confidence"],
    }
    path = os.path.join(args.out, "results", "results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")

    diag = os.path.join(args.out, "results", "l2_diagnostics.json")
    with open(diag, "w", encoding="utf-8") as fh:
        json.dump({k: v for k, v in res2.items()
                   if k not in ("EF_trajectory_eV", "EF_band16_eV", "EF_band84_eV")},
                  fh, indent=2)
        fh.write("\n")

    diag3 = os.path.join(args.out, "results", "l3_diagnostics.json")
    with open(diag3, "w", encoding="utf-8") as fh:
        json.dump({k: v for k, v in res3.items() if k not in ("map", "features")},
                  fh, indent=2)
        fh.write("\n")

    print("\nwrote", path)
    print("wrote", diag)
    print("wrote", diag3)
    for name in ("l1_correction.png", "l2_model_space.png", "l2_trajectory.png",
                 "l2_repeats.png", "l3_map.png", "l3_lcurve.png", "l3_resolution.png"):
        print("wrote", os.path.join(args.out, "figures", name))
    print("C =", round(res["correction_factor"], 4),
          "| n_F =", f"{res['n_fermi_op_cm2']:.3e}", "cm^-2",
          "| mu_F =", round(res["mobility_corrected_cm2_Vs"], 1), "cm^2/Vs")
    print(f"L2: {res2['model']}  Eg = {res2['Eg_kT']:.2f} kT "
          f"(16-84% {eg16:.2f}-{eg84:.2f})  r_gap/r_cone = {res2['ratio']:.3f}")
    print(f"L3: eta = {res3['eta_factor']:.2f} eta0  FWHM = {res3['resolution_probe_pitch']:.2f} pitch  "
          f"delta95 = {res3['detection_limit_fraction']:.2e}  "
          f"a = {res3['a']:.4f}  -> {'PASS' if res3['pass_'] else 'FAIL'} ({res3['confidence']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
