#!/usr/bin/env python3
"""Analyze the periodic buoyancy endpoint correction validation runs."""

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
TESTS = ROOT.parent
NU = 0.1
EOS_RUNS = {
    "thermal": (TESTS / "1DNonlinearEOS/runA_theta", TESTS / "1DNonlinearEOS/runA_theta_periodic_fix"),
    "salty": (TESTS / "1DNonlinearEOS/runB_salt", TESTS / "1DNonlinearEOS/runB_salt_periodic_fix"),
}
PENAL_RUNS = (TESTS / "1DIcePenalization/tau_4dx",
              TESTS / "1DIcePenalization/tau_4dx_periodic_fix_final")
TAU = 2.44140625e-3


def files(directory):
    found = glob.glob(str(directory / "Data_*.h5"))
    return sorted(found, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def solve_bvp(x, rhs, phi=None, tau=None):
    n, dx = len(x), x[1] - x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU / dx**2
    ab[2, :-1] = NU / dx**2
    ab[1, :] = -2.0 * NU / dx**2
    if phi is not None:
        ab[1, :] -= phi / tau
    ab[1, 0] -= NU / dx**2
    ab[1, -1] -= NU / dx**2
    return solve_banded((1, 1), ab, rhs)


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


def scale_error(w, reference):
    scale = float(np.dot(w, reference) / np.dot(reference, reference))
    nominal = float(np.max(np.abs(w-reference)) / np.max(np.abs(w)))
    scaled = float(np.max(np.abs(w-scale*reference)) / np.max(np.abs(w)))
    return scale, nominal, scaled


summary, profiles, history = [], [], []
eos_plot = {}

for name, (old_dir, new_dir) in EOS_RUNS.items():
    old = load(files(old_dir)[-1])
    new_files = files(new_dir)
    new = load(new_files[-1])
    previous = load(new_files[-2])
    b = -np.abs(new["theta"] - 0.4 + 0.5*new["salt"])**2 + 1.75*new["salt"]
    wb = solve_bvp(new["x"], b)
    old_b = -np.abs(old["theta"] - 0.4 + 0.5*old["salt"])**2 + 1.75*old["salt"]
    old_wb = solve_bvp(old["x"], old_b)
    old_scale, old_nominal, _ = scale_error(old["w"], old_wb)
    scale, nominal, scaled = scale_error(new["w"], wb)
    steady = float(np.max(np.abs(new["w"]-previous["w"])) / np.max(np.abs(new["w"])))
    dx = new["x"][1] - new["x"][0]
    b_rec = np.full_like(new["w"], np.nan)
    b_rec[1:-1] = NU*(new["w"][2:]-2*new["w"][1:-1]+new["w"][:-2])/dx**2
    mask = np.abs(b[1:-1]) > 1e-5
    ratio = b_rec[1:-1][mask] / b[1:-1][mask]

    vals = {
        "old_scale": old_scale, "old_nominal_relative_error": old_nominal,
        "fixed_time_final": new["time"], "fixed_scale": scale,
        "fixed_nominal_relative_error": nominal,
        "fixed_scaled_relative_error": scaled,
        "fixed_steady_relative_change": steady,
        "fixed_median_reconstructed_force_ratio": np.median(ratio),
        "fixed_max_u": new["u_max"], "fixed_max_v": new["v_max"],
    }
    if name == "thermal":
        model = lambda th, beta, tmd, q: -beta*np.abs(th-tmd)**q
        fit, _ = curve_fit(model, new["theta"][1:-1], b_rec[1:-1], p0=(1, .4, 2))
        vals.update(fitted_betaT=fit[0], fitted_Tmd0=fit[1], fitted_q=fit[2])
    else:
        c2, c1, c0 = np.polyfit(new["salt"][1:-1], b_rec[1:-1], 2)
        vals.update(fitted_c0=c0, fitted_c1=c1, fitted_c2=c2)
    summary.extend((name, metric, float(value)) for metric, value in vals.items())

    for filename in new_files:
        snap = load(filename)
        bb = -np.abs(snap["theta"]-.4+.5*snap["salt"])**2+1.75*snap["salt"]
        ref = solve_bvp(snap["x"], bb)
        sc, err, _ = scale_error(snap["w"], ref)
        history.append((name, snap["time"], sc, err))
    for row in zip(new["x"], new["w"], wb, b, b_rec):
        profiles.append((name, *map(float, row)))
    eos_plot[name] = (old, old_wb, new, wb)

# Penalization: compare original t=4 result with the refined corrected run.
old_dir, new_dir = PENAL_RUNS
old = load(files(old_dir)[-1])
new_files = files(new_dir)
new = load(new_files[-1])
previous = load(new_files[-2])

def penal_reference(snap):
    phi = np.clip(1.0-snap["F"], 0.0, 1.0)
    return solve_bvp(snap["x"], 2.0*snap["salt"], phi, TAU)

old_ref = penal_reference(old)
new_ref = penal_reference(new)
penal_old, penal_new = old, new
old_scale, old_nominal, _ = scale_error(old["w"], old_ref)
scale, nominal, scaled = scale_error(new["w"], new_ref)
steady = float(np.max(np.abs(new["w"]-previous["w"])) / np.max(np.abs(new["w"])))
dx = new["x"][1]-new["x"][0]
delta = np.sqrt(NU*TAU)
delta_discrete = dx/np.arccosh(1.0+dx**2/(2.0*delta**2))
wa = np.abs(new["w"])
slip = float(abs(np.interp(.5, new["x"], new["w"])))
fit_mask = (new["x"]>.53) & (new["x"]<.8) & (wa>1e-9) & (wa<.8*slip)
delta_fit = float(-1.0/np.polyfit(new["x"][fit_mask], np.log(wa[fit_mask]), 1)[0])
summary.extend(("penalization_4dx", metric, float(value)) for metric, value in {
    "old_scale": old_scale, "old_nominal_relative_error": old_nominal,
    "fixed_time_final": new["time"], "fixed_scale": scale,
    "fixed_nominal_relative_error": nominal,
    "fixed_scaled_relative_error": scaled,
    "fixed_steady_relative_change": steady,
    "delta_discrete": delta_discrete, "delta_fitted": delta_fit,
    "delta_relative_error": abs(delta_fit/delta_discrete-1.0),
    "interface_slip": slip, "fixed_max_u": new["u_max"], "fixed_max_v": new["v_max"],
}.items())
for filename in new_files:
    snap = load(filename)
    ref = penal_reference(snap)
    sc, err, _ = scale_error(snap["w"], ref)
    history.append(("penalization_4dx", snap["time"], sc, err))
for row in zip(new["x"], new["w"], new_ref, new["F"], new["salt"]):
    profiles.append(("penalization_4dx", *map(float, row)))

with open(ROOT/"periodic_fix_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "metric", "value")); w.writerows(summary)
with open(ROOT/"periodic_fix_history.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "time", "scale", "nominal_relative_error")); w.writerows(history)
with open(ROOT/"periodic_fix_profiles.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(("case", "x", "measured_w", "reference_w", "field_1", "field_2")); w.writerows(profiles)

fig, ax = plt.subplots(figsize=(7.2, 4.5))
labels = ["thermal EOS", "salty EOS", "penalization 4dx"]
old_scales = [dict((m,v) for c,m,v in summary if c==key)["old_scale"]
              for key in ("thermal", "salty", "penalization_4dx")]
new_scales = [dict((m,v) for c,m,v in summary if c==key)["fixed_scale"]
              for key in ("thermal", "salty", "penalization_4dx")]
xpos = np.arange(3); width=.35
ax.bar(xpos-width/2, old_scales, width, label="before fix")
ax.bar(xpos+width/2, new_scales, width, label="after fix")
ax.axhline(1, color="k", ls="--"); ax.set_xticks(xpos, labels)
ax.set(ylabel="measured / nominal BVP amplitude", ylim=(0.68, 1.04))
ax.grid(axis="y", alpha=.3); ax.legend(); fig.tight_layout()
fig.savefig(ROOT/"periodic_fix_amplitudes.png", dpi=160); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, name in zip(axes, ("thermal", "salty")):
    old, old_ref, new, ref = eos_plot[name]
    ax.plot(old["x"], old["w"], color="C1", label="before fix")
    ax.plot(new["x"], new["w"], "k-", label="after fix")
    ax.plot(new["x"], ref, "r--", label="nominal discrete BVP")
    ax.set(title=name, xlabel="x", ylabel="w"); ax.grid(alpha=.3); ax.legend()
fig.tight_layout(); fig.savefig(ROOT/"periodic_fix_eos_profiles.png", dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(penal_old["x"], penal_old["w"], color="C1", label="before fix")
ax.plot(penal_new["x"], penal_new["w"], "k-", label="after fix")
ax.plot(penal_new["x"], new_ref, "r--", label="nominal diffuse Brinkman BVP")
ax.axvline(.5, color="gray", ls=":")
ax.set(xlabel="x", ylabel="w", title="4-dx ICE penalization")
ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
fig.savefig(ROOT/"periodic_fix_penalization_profile.png", dpi=160); plt.close(fig)

print("Wrote periodic-fix analysis to", ROOT)
for case, metric, value in summary:
    if metric in ("old_scale", "fixed_scale", "fixed_steady_relative_change",
                  "fixed_median_reconstructed_force_ratio", "delta_relative_error"):
        print(case, metric, f"{value:.12g}")
