#!/usr/bin/env python3
"""Run the frozen pipeline on a blind data set and log what happened.

    python run_blind.py --data <blind dir> --out blind

Writes into <out>:

    blind_results.json   the same schema as results.json, dataset_set = blind
    run_log.txt          one line per phase with wall-clock time and headline
    pipeline_commit.txt  the hash the run belongs to
    run_times.json       machine-readable timings
    figures/             the same seven figures, from the blind data

The pipeline is the same code as reproduce.py; nothing here is specific to the
blind set, which is the only thing worth demonstrating under time pressure.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "src"))
sys.path.insert(0, os.path.join(HERE, ".deps"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mplconfig"))

import l1  # noqa: E402
import l2  # noqa: E402
import l3  # noqa: E402
from common import load_group  # noqa: E402


def pipeline_commit() -> str:
    """Frozen commit of this pipeline, or an explicit marker if there is none."""
    try:
        top = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=HERE,
                             capture_output=True, text=True, timeout=10)
        if top.returncode != 0:
            return "UNCOMMITTED-LOCAL-RUN"
        if os.path.normcase(os.path.realpath(top.stdout.strip())) != \
                os.path.normcase(os.path.realpath(HERE)):
            return "UNCOMMITTED-LOCAL-RUN"
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=HERE,
                              capture_output=True, text=True, timeout=10)
        if head.returncode == 0 and head.stdout.strip():
            return head.stdout.strip()
    except Exception:
        pass
    return "UNCOMMITTED-LOCAL-RUN"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=os.path.join(HERE, "blind"))
    args = ap.parse_args()
    out = args.out
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)

    t_start = time.perf_counter()
    log = []

    def note(line):
        print(line)
        log.append(line)

    note(f"QXU6028 G18 blind run   {time.strftime('%Y-%m-%dT%H:%M:%S')}")
    note(f"python {platform.python_version()}   {platform.platform()}")
    note(f"data   {args.data}")
    note(f"commit {pipeline_commit()}")
    note("")

    data = load_group(args.data)
    params = data["params"]
    note(f"params c = {params['c']:.9f}  T = {params['T']} K  "
         f"E_F,op = {params['EF_op_eV']:.6f} eV  rho = {params['rho_ohm_sq']:.3f} ohm/sq")

    t0 = time.perf_counter()
    res1 = l1.run(params, data["l1"])
    t1 = time.perf_counter()
    note(f"L1  {t1 - t0:6.2f} s   C = {res1['correction_factor']:.4f}  "
         f"n_F = {res1['n_fermi_op_cm2']:.3e} cm^-2  "
         f"mu_F = {res1['mobility_corrected_cm2_Vs']:.1f} cm^2/Vs")

    t0 = time.perf_counter()
    res2 = l2.run(args.data)
    t1 = time.perf_counter()
    eg16, eg84 = res2["Eg_uncertainty_kT"]
    note(f"L2  {t1 - t0:6.2f} s   model = {res2['model']}  Eg = {res2['Eg_kT']:.2f} kT "
         f"(16-84% {eg16:.2f}-{eg84:.2f})  r_gap/r_cone = {res2['ratio']:.3f}  "
         f"gap fraction {res2['gap_fraction']:.3f}")

    t0 = time.perf_counter()
    res3 = l3.run(args.data)
    t1 = time.perf_counter()
    note(f"L3  {t1 - t0:6.2f} s   eta = {res3['eta_factor']:.2f} eta0  "
         f"FWHM = {res3['resolution_probe_pitch']:.2f} pitch  "
         f"delta95 = {res3['detection_limit_fraction']:.2e}  a = {res3['a']:.4f}  "
         f"-> {'PASS' if res3['pass_'] else 'FAIL'} ({res3['confidence']})")

    t0 = time.perf_counter()
    l1.figure(params, os.path.join(figdir, "l1_correction.png"))
    l2.figures(args.data, res2, figdir)
    l3.figures(args.data, res3, figdir)
    t1 = time.perf_counter()
    note(f"fig {t1 - t0:6.2f} s   seven figures written to {figdir}")

    with open(os.path.join(args.data, "blind_results.json"), encoding="utf-8") as fh:
        results = json.load(fh)
    results["pipeline_commit"] = pipeline_commit()
    results["L1"] = {
        "correction_factor": round(res1["correction_factor"], 4),
        "n_fermi_op_cm2": res1["n_fermi_op_cm2"],
        "mobility_corrected_cm2_Vs": round(res1["mobility_corrected_cm2_Vs"], 1),
    }
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

    finite = all(isinstance(v, (int, float)) and v == v and abs(v) < 1e6
                 for v in results["L3"]["map"])
    note(f"checks  map finite: {finite} ({len(results['L3']['map'])} values)  "
         f"group {results['group']}  dataset {results['dataset_set']}")

    path = os.path.join(out, "blind_results.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")
    with open(os.path.join(out, "pipeline_commit.txt"), "w", encoding="utf-8") as fh:
        fh.write(pipeline_commit() + "\n")

    total = time.perf_counter() - t_start
    note(f"total {total:.1f} s")
    if not finite:
        note("STOP: a map value is not finite")
        with open(os.path.join(out, "run_log.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(log) + "\n")
        return 2

    with open(os.path.join(out, "run_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(out, "run_times.json"), "w", encoding="utf-8") as fh:
        json.dump({"total_s": round(total, 1)}, fh, indent=2)
        fh.write("\n")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
