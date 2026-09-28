#!/usr/bin/env python3
"""Audit the fixed-pressure-gradient Poiseuille amplitude diagnostic.

Reads the completed HDF5 outputs from run_exact_all_liquid (256 cells,
initialized with the nominal parabola) and run_zero_all_liquid (64 cells,
initialized at rest).
Writes CSV summaries/profiles/history and PNG plots beside this script.
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
RUNS = ("run_exact_all_liquid", "run_zero_all_liquid")
NU = 0.1
DPDX = -1.0
FORCING = -DPDX


def output_files(run):
    files = glob.glob(str(ROOT / run / "Data_*.h5"))
    return sorted(files, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def discrete_solution(y):
    """Solve nu*u''=dpdx with midpoint no-slip ghost-cell conditions."""
    n = len(y)
    dx = y[1] - y[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU / dx**2
    ab[2, :-1] = NU / dx**2
    ab[1, :] = -2.0 * NU / dx**2
    # Ghost velocity is -u at both walls, giving diagonal coefficient -3.
    ab[1, 0] -= NU / dx**2
    ab[1, -1] -= NU / dx**2
    return solve_banded((1, 1), ab, np.full(n, DPDX))


def load_snapshot(filename):
    with h5py.File(filename, "r") as f:
        y = f["grid/yc"][:-1]
        # The physical y cells exclude the final no-slip ghost.  Periodic
        # copies in x/z are identical; averaging also measures nonuniformity.
        u3 = f["u"][:, :-1, :]
        u = u3.mean(axis=(0, 2))
        v_max = float(np.max(np.abs(f["v"][:])))
        w_max = float(np.max(np.abs(f["w"][:])))
        phase = f["VOF/C_L"][:, :-1, :]
        time = float(np.asarray(f["time"]).reshape(-1)[0])
    reference = discrete_solution(y)
    scale = float(np.dot(u, reference) / np.dot(reference, reference))
    rel_linf = float(np.max(np.abs(u - reference)) / np.max(np.abs(reference)))
    nonuniform = float(np.max(np.abs(u3 - u[None, :, None])))
    return {
        "time": time, "y": y, "u": u, "reference": reference,
        "scale": scale, "rel_linf": rel_linf,
        "v_max": v_max, "w_max": w_max, "nonuniform": nonuniform,
        "phase_min": float(np.min(phase)), "phase_max": float(np.max(phase)),
    }


def reconstructed_force(y, u):
    dx = y[1] - y[0]
    d2 = np.empty_like(u)
    d2[1:-1] = (u[2:] - 2.0*u[1:-1] + u[:-2]) / dx**2
    d2[0] = (u[1] - 3.0*u[0]) / dx**2
    d2[-1] = (u[-2] - 3.0*u[-1]) / dx**2
    return -NU * d2


history = []
profiles = []
summary = []
run_data = {}

for run in RUNS:
    snapshots = [load_snapshot(f) for f in output_files(run)]
    first, previous, last = snapshots[0], snapshots[-2], snapshots[-1]
    force = reconstructed_force(last["y"], last["u"])
    steady_change = np.max(np.abs(last["u"] - previous["u"])) / np.max(np.abs(last["u"]))
    metrics = {
        "time_final": last["time"],
        "grid_cells_y": len(last["y"]),
        "profile_scale_measured_over_discrete_exact": last["scale"],
        "profile_relative_Linf_error": last["rel_linf"],
        "u_max_measured": np.max(last["u"]),
        "u_max_discrete_exact": np.max(last["reference"]),
        "u_bulk_measured": np.mean(last["u"]),
        "u_bulk_discrete_exact": np.mean(last["reference"]),
        "force_reconstructed_mean": np.mean(force),
        "force_reconstructed_median": np.median(force),
        "force_reconstructed_max_abs_error": np.max(np.abs(force - FORCING)),
        "last_output_relative_change": steady_change,
        "max_transverse_v": last["v_max"],
        "max_transverse_w": last["w_max"],
        "streamwise_spanwise_nonuniformity": last["nonuniform"],
        "phase_indicator_min": last["phase_min"],
        "phase_indicator_max": last["phase_max"],
    }
    for metric, value in metrics.items():
        summary.append((run, metric, float(value)))
    for snap in snapshots:
        history.append((run, snap["time"], snap["scale"], snap["rel_linf"],
                        np.max(snap["u"]), np.mean(snap["u"])))
    continuum = FORCING * last["y"] * (1.0-last["y"]) / (2.0*NU)
    for vals in zip(last["y"], last["u"], last["reference"], continuum, force):
        profiles.append((run, *map(float, vals)))
    run_data[run] = (snapshots, force)

# Late-time exponential fit independently estimates the zero-start asymptote.
zero_history = run_data["run_zero_all_liquid"][0]
t = np.array([s["time"] for s in zero_history])
a = np.array([s["scale"] for s in zero_history])
late = t >= 2.0
model = lambda tt, asymptote, amplitude, decay: asymptote - amplitude*np.exp(-decay*tt)
fit, _ = curve_fit(model, t[late], a[late], p0=(1.0, 1.0, np.pi**2*NU), maxfev=10000)
summary.extend([
    ("run_zero_all_liquid_fit", "fitted_asymptotic_scale", float(fit[0])),
    ("run_zero_all_liquid_fit", "fitted_decay_rate", float(fit[2])),
    ("run_zero_all_liquid_fit", "theoretical_slowest_decay_rate", float(np.pi**2*NU)),
    ("run_zero_all_liquid_fit", "decay_rate_relative_error", float(abs(fit[2]/(np.pi**2*NU)-1.0))),
])

with open(ROOT / "poiseuille_results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(("run", "metric", "value"))
    writer.writerows(summary)

with open(ROOT / "poiseuille_history.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(("run", "time", "profile_scale", "relative_Linf_error", "u_max", "u_bulk"))
    writer.writerows(history)

with open(ROOT / "poiseuille_profiles.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(("run", "y", "u_measured", "u_discrete_exact", "u_continuum", "force_reconstructed"))
    writer.writerows(profiles)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, run in zip(axes, RUNS):
    snapshots, _ = run_data[run]
    last = snapshots[-1]
    continuum = FORCING * last["y"] * (1.0-last["y"]) / (2.0*NU)
    ax.plot(last["y"], last["u"], "k-", lw=2, label="measured")
    ax.plot(last["y"], last["reference"], "r--", label="discrete exact")
    ax.plot(last["y"], continuum, "b:", label="continuum")
    ax.set(title=f"{run}, scale={last['scale']:.8f}", xlabel="y", ylabel="u")
    ax.grid(alpha=.3)
    ax.legend()
fig.tight_layout()
fig.savefig(ROOT / "poiseuille_profiles.png", dpi=160)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
for run in RUNS:
    snapshots, _ = run_data[run]
    ax.plot([s["time"] for s in snapshots], [s["scale"] for s in snapshots], "o-", label=run)
ax.axhline(1.0, color="k", ls="--", label="nominal amplitude = 1")
ax.axhline(0.75, color="r", ls=":", label="EOS/penalization amplitude = 0.75")
tt = np.linspace(2.0, max(t), 200)
ax.plot(tt, model(tt, *fit), color="C1", ls="--", label=f"zero-start fit: asymptote={fit[0]:.9f}")
ax.set(xlabel="time", ylabel="measured / exact profile amplitude", ylim=(-.03, 1.05))
ax.grid(alpha=.3)
ax.legend()
fig.tight_layout()
fig.savefig(ROOT / "poiseuille_amplitude_history.png", dpi=160)
plt.close(fig)

print("Poiseuille audit written to", ROOT)
for run, metric, value in summary:
    if metric in ("profile_scale_measured_over_discrete_exact", "profile_relative_Linf_error",
                  "force_reconstructed_mean", "fitted_asymptotic_scale", "fitted_decay_rate"):
        print(run, metric, f"{value:.12g}")
