#!/usr/bin/env python3
"""Independent completed-run audit for the EOS_NONLINEAR gate.

Reads the existing HDF5 snapshots and writes audit_eos_results.csv,
audit_eos_profiles.csv, audit_eos_buoyancy.png, and audit_eos_velocity.png.
It does not modify simulation inputs, source code, or the original analyzer.
"""

from pathlib import Path
import csv
import glob
import re

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parent
RUNS = ("runA_theta", "runB_salt")
RE, NU = 10.0, 0.1
BETAT, BETAS, TMD0, SLOPE, Q = 1.0, 1.75, 0.4, -0.5, 2.0


def outputs(run):
    files = glob.glob(str(ROOT / run / "Data_*.h5"))
    return sorted(files, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def load(fn):
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        x = f["grid/xc"][:-1]
        th = f["Conc/0"][1:-1, 1:-1, :-1].mean((0, 1))
        s = f["Conc/1"][1:-1, 1:-1, :-1].mean((0, 1))
        F = f["VOF/C_L"][1:-1, 1:-1, :-1].mean((0, 1))
        w3 = f["w"][1:-1, 1:-1, :-1]
        w = w3.mean((0, 1))
        u_max = float(np.max(np.abs(f["u"][:])))
        v_max = float(np.max(np.abs(f["v"][:])))
    return {"t": t, "x": x, "th": th, "s": s, "F": F, "w": w,
            "w_nonuniform": float(np.max(np.abs(w3 - w[None, None, :]))),
            "u_max": u_max, "v_max": v_max}


def eos(th, s):
    return -BETAT * np.abs(th - TMD0 - SLOPE * s) ** Q + BETAS * s


def solve_bvp(x, rhs):
    """Solve nu*w''=rhs with midpoint no-slip walls."""
    n, dx = len(x), x[1] - x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU / dx**2
    ab[2, :-1] = NU / dx**2
    ab[1, :] = -2.0 * NU / dx**2
    ab[1, 0] -= NU / dx**2
    ab[1, -1] -= NU / dx**2
    return solve_banded((1, 1), ab, rhs)


summary, profiles, data = [], [], {}
for run in RUNS:
    fs = outputs(run)
    first, prev, last = load(fs[0]), load(fs[-2]), load(fs[-1])
    x, w = last["x"], last["w"]
    dx = x[1] - x[0]
    b = eos(last["th"], last["s"])
    wb = solve_bvp(x, b)
    scale = float(np.dot(w, wb) / np.dot(wb, wb))
    scaled_error = float(np.max(np.abs(w - scale * wb)) / np.max(np.abs(w)))
    nominal_error = float(np.max(np.abs(w - wb)) / np.max(np.abs(w)))
    b_rec = np.full_like(w, np.nan)
    b_rec[1:-1] = NU * (w[2:] - 2 * w[1:-1] + w[:-2]) / dx**2
    ratio = b_rec[1:-1] / b[1:-1]
    finite = np.isfinite(ratio) & (np.abs(b[1:-1]) > 1e-5)

    values = {
        "t_final": last["t"],
        "steady_relative_change": np.max(np.abs(w - prev["w"])) / np.max(np.abs(w)),
        "F_drift": np.max(np.abs(last["F"] - first["F"])),
        "theta_drift": np.max(np.abs(last["th"] - first["th"])),
        "salinity_drift": np.max(np.abs(last["s"] - first["s"])),
        "bvp_scale_measured_over_nominal": scale,
        "nominal_bvp_relative_error": nominal_error,
        "scaled_bvp_relative_error": scaled_error,
        "median_b_reconstructed_over_input": np.median(ratio[finite]),
        "max_u": last["u_max"], "max_v": last["v_max"],
        "w_nonuniformity": last["w_nonuniform"],
    }

    if run == "runA_theta":
        model = lambda th, beta, tmd, q: -beta * np.abs(th - tmd) ** q
        fit, _ = curve_fit(model, last["th"][1:-1], b_rec[1:-1],
                           p0=(0.75, 0.4, 2.0))
        values.update(fitted_betaT=fit[0], fitted_Tmd0=fit[1], fitted_q=fit[2])
    else:
        c2, c1, c0 = np.polyfit(last["s"][1:-1], b_rec[1:-1], 2)
        values.update(fitted_c0=c0, fitted_c1=c1, fitted_c2=c2)

    for metric, value in values.items():
        summary.append((run, metric, float(value)))
    for i in range(len(x)):
        profiles.append((run, x[i], last["th"][i], last["s"][i], last["F"][i],
                         w[i], wb[i], scale * wb[i], b[i], b_rec[i]))
    data[run] = (last, b, b_rec, wb, scale)

# Matching 1-rank/4-rank comparison at t=1.
one = ROOT / "runA_theta" / "Data_2.h5"
four = ROOT / "runA_theta_4rank" / "Data_2.h5"
if one.exists() and four.exists():
    with h5py.File(one, "r") as a, h5py.File(four, "r") as b:
        for key in ("Conc/0", "Conc/1", "VOF/C_L", "u", "v", "w"):
            summary.append(("mpi_1_vs_4", key, float(np.max(np.abs(a[key][:] - b[key][:])))))

with open(ROOT / "audit_eos_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("run", "metric", "value")); w.writerows(summary)
with open(ROOT / "audit_eos_profiles.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(("run", "x", "theta", "salinity", "F", "w_measured",
                "w_nominal_bvp", "w_scaled_bvp", "b_input", "b_reconstructed"))
    w.writerows(profiles)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, run in zip(axes, RUNS):
    last, b, b_rec, _, scale = data[run]
    ax.plot(last["x"], b, "r-", label="specified EOS b")
    ax.plot(last["x"][1:-1], b_rec[1:-1], "k.", ms=3, label="reconstructed from w")
    ax.plot(last["x"], scale * b, "b--", label=f"{scale:.6f} x specified b")
    ax.set(title=run, xlabel="x", ylabel="buoyancy"); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(ROOT / "audit_eos_buoyancy.png", dpi=160); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, run in zip(axes, RUNS):
    last, _, _, wb, scale = data[run]
    ax.plot(last["x"], last["w"], "k-", label="measured")
    ax.plot(last["x"], wb, "r--", label="nominal BVP")
    ax.plot(last["x"], scale * wb, "b:", lw=2, label=f"scaled BVP ({scale:.6f})")
    ax.set(title=run, xlabel="x", ylabel="w"); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(ROOT / "audit_eos_velocity.png", dpi=160); plt.close(fig)

print("EOS audit written to", ROOT)
for row in summary:
    if row[1] in ("bvp_scale_measured_over_nominal", "scaled_bvp_relative_error",
                  "median_b_reconstructed_over_input"):
        print(row)
