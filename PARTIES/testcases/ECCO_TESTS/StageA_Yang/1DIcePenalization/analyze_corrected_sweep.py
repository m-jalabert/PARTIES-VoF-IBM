#!/usr/bin/env python3
"""Analyze the complete corrected 2/4/8-dx ICE-penalization sweep."""

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
NU = 0.1
X_INT = 0.5
CN = 0.0029296875
RUNS = {
    "tau_2dx": ("tau_2dx_periodic_fix_final", 6.103515625e-4),
    "tau_4dx": ("tau_4dx_periodic_fix_final", 2.44140625e-3),
    "tau_8dx": ("tau_8dx_periodic_fix_final", 9.765625e-3),
}


def outputs(directory):
    found = glob.glob(str(ROOT / directory / "Data_*.h5"))
    return sorted(found, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def load(filename):
    with h5py.File(filename, "r") as h5:
        return {
            "time": float(h5["time"][0]),
            "x": h5["grid/xc"][:-1],
            "theta": h5["Conc/0"][1:-1, 1:-1, :-1].mean((0, 1)),
            "salt": h5["Conc/1"][1:-1, 1:-1, :-1].mean((0, 1)),
            "F": h5["VOF/C_L"][1:-1, 1:-1, :-1].mean((0, 1)),
            "w": h5["w"][1:-1, 1:-1, :-1].mean((0, 1)),
            "u_max": float(np.max(np.abs(h5["u"][:]))),
            "v_max": float(np.max(np.abs(h5["v"][:]))),
        }


def solve_bvp(x, phi, tau, buoyancy):
    n, dx = len(x), x[1]-x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU/dx**2
    ab[2, :-1] = NU/dx**2
    ab[1, :] = -2.0*NU/dx**2 - phi/tau
    ab[1, 0] -= NU/dx**2
    ab[1, -1] -= NU/dx**2
    return solve_banded((1, 1), ab, buoyancy)


def reference(snapshot, tau):
    return solve_bvp(snapshot["x"], np.clip(1.0-snapshot["F"], 0.0, 1.0),
                     tau, 2.0*snapshot["salt"])


def comparison(w, wb):
    scale = float(np.dot(w, wb)/np.dot(wb, wb))
    nominal = float(np.max(np.abs(w-wb))/np.max(np.abs(w)))
    scaled = float(np.max(np.abs(w-scale*wb))/np.max(np.abs(w)))
    return scale, nominal, scaled


summary, profiles, history, plot_data = [], [], [], {}
taus, slips, bvp_slips, sharp_slips = [], [], [], []

for name, (directory, tau) in RUNS.items():
    fs = outputs(directory)
    # Always use the original t=4 gate state for the before-fix comparison;
    # the final 4-dx continuation directory starts from a later refined restart.
    first = load(ROOT/name/"Data_20.h5")
    previous, last = load(fs[-2]), load(fs[-1])
    wb0, wb = reference(first, tau), reference(last, tau)
    old_scale, _, _ = comparison(first["w"], wb0)
    scale, nominal, scaled = comparison(last["w"], wb)
    steady = float(np.max(np.abs(last["w"]-previous["w"]))/np.max(np.abs(last["w"])))
    x, w = last["x"], last["w"]
    dx = x[1]-x[0]
    delta = np.sqrt(NU*tau)
    delta_discrete = dx/np.arccosh(1.0+dx**2/(2.0*delta**2))
    wa = np.abs(w)
    slip = float(abs(np.interp(X_INT, x, w)))
    bvp_slip = float(abs(np.interp(X_INT, x, wb)))
    fit_mask = (x>X_INT+.03) & (x<.8) & (wa>1e-9) & (wa<.8*slip)
    delta_fit = float(-1.0/np.polyfit(x[fit_mask], np.log(wa[fit_mask]), 1)[0])
    band = 6.0*2.0*np.sqrt(2.0)*CN
    far = x>X_INT+band+5.0*delta
    leakage = float(np.max(wa[far])/np.max(wa)) if np.any(far) else 0.0
    X = X_INT
    sharp = abs(X*((2.0*X/NU)*(2.0*delta+X)/(2.0*(delta+X))-2.0*X/(2.0*NU)))

    passed = (abs(scale-1.0)<.01 and nominal<.01 and scaled<1e-3 and
              steady<1e-3 and abs(delta_fit/delta_discrete-1.0)<.02 and leakage<.01)
    values = {
        "old_scale": old_scale, "time_final": last["time"],
        "fixed_scale": scale, "fixed_scale_error": abs(scale-1.0),
        "nominal_bvp_relative_error": nominal,
        "scaled_bvp_relative_error": scaled,
        "steady_relative_change": steady,
        "delta_continuum": delta, "delta_discrete": delta_discrete,
        "delta_fitted": delta_fit,
        "delta_error_vs_discrete": abs(delta_fit/delta_discrete-1.0),
        "delta_error_vs_continuum": abs(delta_fit/delta-1.0),
        "interface_slip_measured": slip, "interface_slip_bvp": bvp_slip,
        "far_ice_leakage": leakage,
        "F_change_from_t4": np.max(np.abs(last["F"]-first["F"])),
        "salt_change_from_t4": np.max(np.abs(last["salt"]-first["salt"])),
        "max_u": last["u_max"], "max_v": last["v_max"],
        "case_pass": int(passed),
    }
    summary.extend((name, metric, float(value)) for metric, value in values.items())
    for filename in fs:
        snap = load(filename)
        ref = reference(snap, tau)
        sc, err, _ = comparison(snap["w"], ref)
        history.append((name, snap["time"], sc, err))
    for row in zip(x, last["F"], last["salt"], w, wb):
        profiles.append((name, tau, *map(float, row)))
    plot_data[name] = (last, wb, delta)
    taus.append(tau); slips.append(slip); bvp_slips.append(bvp_slip); sharp_slips.append(sharp)

taus = np.asarray(taus)
slips = np.asarray(slips)
bvp_slips = np.asarray(bvp_slips)
sharp_slips = np.asarray(sharp_slips)
slope_measured = float(np.polyfit(np.log(taus), np.log(slips), 1)[0])
slope_bvp = float(np.polyfit(np.log(taus), np.log(bvp_slips), 1)[0])
slope_sharp = float(np.polyfit(np.log(taus), np.log(sharp_slips), 1)[0])
sweep_pass = abs(slope_measured-slope_bvp)<.02
summary.extend([
    ("sweep", "slip_exponent_measured", slope_measured),
    ("sweep", "slip_exponent_diffuse_bvp", slope_bvp),
    ("sweep", "slip_exponent_sharp_reference", slope_sharp),
    ("sweep", "slip_exponent_error_vs_diffuse_bvp", abs(slope_measured-slope_bvp)),
    ("sweep", "sweep_pass", float(sweep_pass)),
])

# Corrected-source rank-count check from the same t=4 restart to t=4.2.
mpi1 = ROOT/"tau_4dx_periodic_fix_mpi1/Data_21.h5"
mpi4 = ROOT/"tau_4dx_periodic_fix_mpi4/Data_21.h5"
mpi_worst = 0.0
with h5py.File(mpi1, "r") as one, h5py.File(mpi4, "r") as four:
    for field in ("Conc/0", "Conc/1", "VOF/C_L", "u", "v", "w", "p"):
        difference = float(np.max(np.abs(one[field][:]-four[field][:])))
        mpi_worst = max(mpi_worst, difference)
        summary.append(("mpi_1_vs_4", field, difference))
summary.append(("mpi_1_vs_4", "worst_field_difference", mpi_worst))
summary.append(("mpi_1_vs_4", "mpi_pass", float(mpi_worst < 1e-10)))

with open(ROOT/"corrected_penalization_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "metric", "value")); w.writerows(summary)
with open(ROOT/"corrected_penalization_profiles.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "darcy_tau", "x", "F", "salt", "w_measured", "w_bvp")); w.writerows(profiles)
with open(ROOT/"corrected_penalization_history.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "time", "scale", "nominal_relative_error")); w.writerows(history)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
for name in RUNS:
    last, wb, delta = plot_data[name]
    axes[0].plot(last["x"], last["w"], label=f"{name} measured")
    axes[0].plot(last["x"], wb, "k:", lw=1)
    axes[1].semilogy(last["x"], np.abs(last["w"])+1e-16, label=f"{name}, delta={delta:.5f}")
for ax in axes:
    ax.axvline(X_INT, color="gray", ls="--"); ax.grid(alpha=.3); ax.legend()
axes[0].set(xlabel="x", ylabel="w", title="Corrected profiles; dotted = nominal BVP")
axes[1].set(xlabel="x", ylabel="|w|", title="Brinkman decay in ice", ylim=(1e-12, 1))
fig.tight_layout(); fig.savefig(ROOT/"corrected_penalization_profiles.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 4.6))
ax.loglog(taus, slips, "ko-", label=f"measured, slope={slope_measured:.3f}")
ax.loglog(taus, bvp_slips, "bs--", label=f"diffuse BVP, slope={slope_bvp:.3f}")
ax.loglog(taus, sharp_slips, "r^-.", label=f"sharp reference, slope={slope_sharp:.3f}")
ax.set(xlabel="darcy_tau", ylabel="|w(X)|", title="Corrected interface-slip scaling")
ax.grid(alpha=.3, which="both"); ax.legend(); fig.tight_layout()
fig.savefig(ROOT/"corrected_penalization_slip_scaling.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.6))
for name in RUNS:
    rows = [(t, sc) for case, t, sc, _ in history if case==name]
    ax.plot([r[0] for r in rows], [r[1] for r in rows], "o-", label=name)
ax.axhline(1.0, color="k", ls="--")
ax.set(xlabel="time", ylabel="measured / nominal BVP amplitude", title="Relaxation after periodic forcing fix")
ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
fig.savefig(ROOT/"corrected_penalization_history.png", dpi=160); plt.close(fig)

print("Corrected penalization sweep written to", ROOT)
for row in summary:
    if row[1] in ("fixed_scale", "steady_relative_change", "delta_error_vs_discrete",
                  "case_pass", "slip_exponent_measured", "slip_exponent_diffuse_bvp",
                  "sweep_pass", "worst_field_difference", "mpi_pass"):
        print(row)
