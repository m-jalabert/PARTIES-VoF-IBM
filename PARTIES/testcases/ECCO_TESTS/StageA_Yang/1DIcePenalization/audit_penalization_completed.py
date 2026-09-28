#!/usr/bin/env python3
"""Independent completed-run audit for ICE_PENALIZATION.

Uses the correct -z gravity sign and the measured diffuse F/s fields. Writes
audit_penalization_results.csv, audit_penalization_profiles.csv, and three PNG
figures without changing inputs, solver code, or the original analyzer.
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

ROOT = Path(__file__).resolve().parent
RE, NU, X_INT, CN = 10.0, 0.1, 0.5, 0.0029296875
RUNS = {"tau_2dx": 6.103515625e-4,
        "tau_4dx": 2.44140625e-3,
        "tau_8dx": 9.765625e-3}


def outputs(run):
    fs = glob.glob(str(ROOT / run / "Data_*.h5"))
    return sorted(fs, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def load(fn):
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        x = f["grid/xc"][:-1]
        th = f["Conc/0"][1:-1, 1:-1, :-1].mean((0, 1))
        s = f["Conc/1"][1:-1, 1:-1, :-1].mean((0, 1))
        F = f["VOF/C_L"][1:-1, 1:-1, :-1].mean((0, 1))
        w = f["w"][1:-1, 1:-1, :-1].mean((0, 1))
        u_max = float(np.max(np.abs(f["u"][:])))
        v_max = float(np.max(np.abs(f["v"][:])))
    return {"t": t, "x": x, "th": th, "s": s, "F": F, "w": w,
            "u_max": u_max, "v_max": v_max}


def solve_bvp(x, phi, tau, b):
    """Solve nu*w''-(phi/tau)w=b for gravity directed along -z."""
    n, dx = len(x), x[1] - x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU / dx**2
    ab[2, :-1] = NU / dx**2
    ab[1, :] = -2.0 * NU / dx**2 - phi / tau
    ab[1, 0] -= NU / dx**2
    ab[1, -1] -= NU / dx**2
    return solve_banded((1, 1), ab, b)


summary, profiles, data = [], [], {}
slips, bvp_slips, sharp_slips, taus = [], [], [], []
for run, tau in RUNS.items():
    fs = outputs(run)
    first, prev, last = load(fs[0]), load(fs[-2]), load(fs[-1])
    x, w = last["x"], last["w"]
    dx = x[1] - x[0]
    phi = np.clip(1.0 - last["F"], 0.0, 1.0)
    b = 2.0 * last["s"]
    wb = solve_bvp(x, phi, tau, b)
    scale = float(np.dot(w, wb) / np.dot(wb, wb))
    scaled_error = float(np.max(np.abs(w - scale * wb)) / np.max(np.abs(w)))
    nominal_error = float(np.max(np.abs(w - wb)) / np.max(np.abs(w)))
    delta = np.sqrt(NU * tau)
    delta_disc = dx / np.arccosh(1.0 + dx**2 / (2.0 * delta**2))

    # Fit before the numerical floor/east-wall reflection dominates.
    wa = np.abs(w); slip = float(abs(np.interp(X_INT, x, w)))
    fit_mask = (x > X_INT + 0.03) & (x < 0.8) & (wa > 1e-9) & (wa < 0.8 * slip)
    decay = float(-1.0 / np.polyfit(x[fit_mask], np.log(wa[fit_mask]), 1)[0])
    band = 6.0 * 2.0 * np.sqrt(2.0) * CN
    far = x > X_INT + band + 5.0 * delta
    leakage = float(np.max(wa[far]) / np.max(wa)) if np.any(far) else 0.0

    X = X_INT
    sharp = abs(X * ((2.0 * X / NU) * (2 * delta + X) /
                         (2 * (delta + X)) - 2.0 * X / (2 * NU)))
    bvp_slip = float(abs(np.interp(X_INT, x, wb)))
    slips.append(slip); bvp_slips.append(bvp_slip); sharp_slips.append(sharp); taus.append(tau)

    values = {
        "t_final": last["t"],
        "steady_relative_change": np.max(np.abs(w - prev["w"])) / np.max(np.abs(w)),
        "F_drift": np.max(np.abs(last["F"] - first["F"])),
        "theta_drift": np.max(np.abs(last["th"] - first["th"])),
        "salinity_drift": np.max(np.abs(last["s"] - first["s"])),
        "bvp_scale_measured_over_nominal": scale,
        "nominal_bvp_relative_error": nominal_error,
        "scaled_bvp_relative_error": scaled_error,
        "delta_continuum": delta, "delta_discrete": delta_disc,
        "delta_fitted": decay,
        "delta_error_vs_discrete": abs(decay - delta_disc) / delta_disc,
        "delta_error_vs_continuum": abs(decay - delta) / delta,
        "interface_slip_measured": slip, "interface_slip_bvp": bvp_slip,
        "interface_slip_sharp": sharp, "far_ice_leakage": leakage,
        "max_u": last["u_max"], "max_v": last["v_max"],
    }
    for metric, value in values.items(): summary.append((run, metric, float(value)))
    for i in range(len(x)):
        profiles.append((run, tau, x[i], last["F"][i], last["s"][i], w[i],
                         wb[i], scale * wb[i]))
    data[run] = (last, wb, scale, delta)

taus = np.asarray(taus); slips = np.asarray(slips); bvp_slips = np.asarray(bvp_slips)
sharp_slips = np.asarray(sharp_slips)
slope_measured = float(np.polyfit(np.log(taus), np.log(slips), 1)[0])
slope_bvp = float(np.polyfit(np.log(taus), np.log(bvp_slips), 1)[0])
slope_sharp = float(np.polyfit(np.log(taus), np.log(sharp_slips), 1)[0])
summary += [("sweep", "slip_exponent_measured", slope_measured),
            ("sweep", "slip_exponent_diffuse_bvp", slope_bvp),
            ("sweep", "slip_exponent_sharp_reference", slope_sharp)]

one = ROOT / "tau_4dx" / "Data_2.h5"
four = ROOT / "tau_4dx_4rank" / "Data_2.h5"
if one.exists() and four.exists():
    with h5py.File(one, "r") as a, h5py.File(four, "r") as b:
        for key in ("Conc/0", "Conc/1", "VOF/C_L", "u", "v", "w"):
            summary.append(("mpi_1_vs_4", key, float(np.max(np.abs(a[key][:] - b[key][:])))))

with open(ROOT / "audit_penalization_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("run", "metric", "value")); w.writerows(summary)
with open(ROOT / "audit_penalization_profiles.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(("run", "darcy_tau", "x", "F", "salinity", "w_measured",
                "w_nominal_correct_sign_bvp", "w_scaled_bvp"))
    w.writerows(profiles)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for run in RUNS:
    last, wb, scale, delta = data[run]
    axes[0].plot(last["x"], last["w"], label=f"{run} measured")
    axes[0].plot(last["x"], scale * wb, "k:", lw=1)
    axes[1].semilogy(last["x"], np.abs(last["w"]) + 1e-16,
                    label=f"{run}, delta={delta:.5f}")
axes[0].axvline(X_INT, color="gray", ls="--"); axes[1].axvline(X_INT, color="gray", ls="--")
axes[0].set(xlabel="x", ylabel="w", title="Measured profiles; dotted = scaled diffuse BVP")
axes[1].set(xlabel="x", ylabel="|w|", title="Brinkman decay in ice", ylim=(1e-12, 1))
for ax in axes: ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(ROOT / "audit_penalization_profiles.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.6))
for run in RUNS:
    last, _, _, delta = data[run]
    distance = last["x"] - X_INT
    ax.semilogy(distance, np.abs(last["w"]) + 1e-16, label=f"{run}")
ax.set(xlabel="x - X_interface", ylabel="|w|", title="Measured exponential decay lengths")
ax.set_xlim(0, .3); ax.set_ylim(1e-12, .2); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(ROOT / "audit_penalization_decay.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.6))
ax.loglog(taus, slips, "ko-", label=f"measured, slope={slope_measured:.3f}")
ax.loglog(taus, bvp_slips, "bs--", label=f"diffuse BVP, slope={slope_bvp:.3f}")
ax.loglog(taus, sharp_slips, "r^-.", label=f"sharp reference, slope={slope_sharp:.3f}")
ax.set(xlabel="darcy_tau", ylabel="|w(X)|", title="Interface-slip scaling")
ax.grid(alpha=.3, which="both"); ax.legend()
fig.tight_layout(); fig.savefig(ROOT / "audit_penalization_slip_scaling.png", dpi=160); plt.close(fig)

print("Penalization audit written to", ROOT)
for row in summary:
    if row[1] in ("bvp_scale_measured_over_nominal", "scaled_bvp_relative_error",
                  "delta_fitted", "slip_exponent_measured", "slip_exponent_diffuse_bvp"):
        print(row)
