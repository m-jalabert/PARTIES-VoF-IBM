#!/usr/bin/env python3
"""
EOS_NONLINEAR validation gate - roadmap section A.3c.

Steady unidirectional buoyancy-driven flow in an all-liquid box: frozen
scalar profiles theta(x), s(x) (Pe = 1e12), gravity along the thin periodic
z direction, x = no-slip walls, y = free-slip walls.  The steady flow w(x)
satisfies the DISCRETE balance exactly (u = v = 0, dp = 0, advection = 0):

    (1/Re) * (w[i+1] - 2 w[i] + w[i-1]) / dx^2  =  b(theta_i, s_i)

so the discrete curvature of the measured w(x) reconstructs the Roquet
buoyancy pointwise:

    b(theta, s) = -betaT * |theta - Tmd0 - slope*s|^q + betaS * s

Sub-runs (same [eos] block: q=2, betaT=1, betaS=1.75, Tmd0=0.4, slope=-0.5):
  runA_theta: theta = erf ramp 0->1, s = 0    -> b = -(theta-0.4)^2
              certifies betaT, q, Tmd0
  runB_salt : theta = 1, s = erf ramp 0->1    -> b = -(0.6+0.5 s)^2 + 1.75 s
              certifies eos_Tmd_slope (magnitude AND sign) and betaS

Gate criteria (all must hold in BOTH sub-runs, exit 0 = PASS):
  1. steadiness:   max|w(t_end) - w(t_end - 0.1)| <= 1e-5 * max|w|
  2. frozen state: F = 1 exactly (<=1e-12); theta, s drift <= 1e-9
  3. profile:      max|w - w_BVP| <= 1e-6  (same discrete operator, solved
                   directly as a boundary-value problem from the inputs)
  4. pointwise EOS: max|b_rec - b_input| <= 1e-5 over interior cells
  5. fits:  runA   betaT (1%), Tmd0 (+-0.005), q (+-0.02) from a 3-parameter
                   fit of b_rec(theta);
            runB   quadratic fit b_rec(s) = c0 + c1 s + c2 s^2 against
                   (-0.36, +1.15, -0.25): c1 flips to +2.35 if the SIGN of
                   eos_Tmd_slope is miswired -> sign-discriminating
  6. transverse quiescence: max|u|, max|v| <= 1e-12
  7. [optional, if runA_theta_4rank exists] MPI: 1-rank vs 4-rank fields at
     the matching output agree to <= 1e-13

Writes eos_results.csv (gate summary), eos_profiles.csv, fig_*.png.
Run from testcases/1DNonlinearEOS/ after both sub-runs: python3 analyze_eos.py
"""

import csv
import glob
import os
import re
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.linalg import solve_banded
from scipy.optimize import curve_fit

# --- parameters (must match parties.inp) --------------------------------------
RE     = 10.0
NU     = 1.0 / RE
OUT_DT = 0.5              # output_time_interval (for the MPI-check index)
Q      = 2.0
BETAT  = 1.0
BETAS  = 1.75
TMD0   = 0.4
SLOPE  = -0.5

TOL_STEADY   = 1.0e-5     # relative, successive outputs
TOL_FROZEN_F = 1.0e-12
TOL_FROZEN_C = 1.0e-9
TOL_BVP      = 1.0e-6     # absolute on w
TOL_POINTWISE= 1.0e-5     # absolute on b
TOL_BETAT    = 0.01       # relative
TOL_TMD0     = 0.005      # absolute
TOL_Q        = 0.02       # absolute
TOL_C_RUNB   = 0.01       # absolute on quadratic-fit coefficients
TOL_UV       = 1.0e-12
TOL_MPI      = 1.0e-13

def b_eos(theta, s):
    dth = theta - TMD0 - SLOPE * s
    return -BETAT * np.abs(dth) ** Q + BETAS * s


def load_snapshot(fn):
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        xc = f["grid/xc"][:]
        nx = len(xc) - 1
        th = f["Conc/0"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        s  = f["Conc/1"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        F  = f["VOF/C_L"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        w  = f["w"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        wdev = np.abs(f["w"][:][1:-1, 1:-1, :nx] - w[None, None, :]).max()
        umax = np.abs(f["u"][:]).max()
        vmax = np.abs(f["v"][:]).max()
    return dict(t=t, x=xc[:nx], th=th, s=s, F=F, w=w, wdev=wdev,
                umax=umax, vmax=vmax)


def sorted_outputs(d):
    files = sorted(glob.glob(os.path.join(d, "Data_*.h5")),
                   key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))
    if len(files) < 3:
        sys.exit(f"{d}: only {len(files)} outputs - run the case first")
    return files


def solve_bvp(x, rhs):
    """Solve nu*w'' = rhs with w = 0 on both walls (half-cell ghost:
    w_ghost = -w_interior), same stencil as the code's uniform-mu operator."""
    n = len(x)
    dx = x[1] - x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = 1.0                # upper diag
    ab[2, :-1] = 1.0               # lower diag
    ab[1, :] = -2.0
    ab[1, 0] = -3.0                # ghost = -w0  (wall midpoint zero)
    ab[1, -1] = -3.0
    return solve_banded((1, 1), ab, rhs * dx * dx / NU)


checks = []   # (run, criterion, value, threshold, pass)
rows_profiles = []
fit_results = {}

runs = ["runA_theta", "runB_salt"]
snap = {}

for run in runs:
    files = sorted_outputs(run)
    first = load_snapshot(files[0])
    prev = load_snapshot(files[-2])
    last = load_snapshot(files[-1])
    snap[run] = (first, prev, last)

    x, dx = last["x"], last["x"][1] - last["x"][0]
    wmax = np.abs(last["w"]).max()

    # 1. steadiness
    dsteady = np.abs(last["w"] - prev["w"]).max() / wmax
    checks.append((run, "steadiness d|w|/max|w| (last two outputs)",
                   dsteady, TOL_STEADY, dsteady <= TOL_STEADY))

    # 2. frozen background state
    dF = np.abs(last["F"] - 1.0).max()
    dth = np.abs(last["th"] - first["th"]).max()
    ds = np.abs(last["s"] - first["s"]).max()
    checks.append((run, "all-liquid F=1 max|F-1|", dF, TOL_FROZEN_F, dF <= TOL_FROZEN_F))
    checks.append((run, "frozen theta drift", dth, TOL_FROZEN_C, dth <= TOL_FROZEN_C))
    checks.append((run, "frozen s drift", ds, TOL_FROZEN_C, ds <= TOL_FROZEN_C))

    # 3. w(x) against the direct BVP solution of the same discrete operator
    b_in = b_eos(last["th"], last["s"])
    w_bvp = solve_bvp(x, b_in)
    dbvp = np.abs(last["w"] - w_bvp).max()
    checks.append((run, "max|w - w_BVP|", dbvp, TOL_BVP, dbvp <= TOL_BVP))

    # 4. pointwise buoyancy reconstruction (interior cells)
    w = last["w"]
    b_rec = np.full_like(w, np.nan)
    b_rec[1:-1] = NU * (w[2:] - 2.0 * w[1:-1] + w[:-2]) / dx**2
    err_b = np.abs(b_rec[1:-1] - b_in[1:-1]).max()
    checks.append((run, "pointwise max|b_rec - b_input|",
                   err_b, TOL_POINTWISE, err_b <= TOL_POINTWISE))

    # 6. transverse quiescence + z-uniformity of w
    checks.append((run, "max|u|", last["umax"], TOL_UV, last["umax"] <= TOL_UV))
    checks.append((run, "max|v|", last["vmax"], TOL_UV, last["vmax"] <= TOL_UV))
    checks.append((run, "w z/y-nonuniformity", last["wdev"], TOL_UV,
                   last["wdev"] <= TOL_UV))

    # 5. physics fits
    if run == "runA_theta":
        def model(theta, bT, tmd, q):
            return -bT * np.abs(theta - tmd) ** q
        m = slice(1, -1)
        (bT_f, tmd_f, q_f), _ = curve_fit(model, last["th"][m], b_rec[m],
                                          p0=[0.8, 0.3, 1.8])
        checks.append((run, "fit betaT", bT_f, f"{BETAT}+-{TOL_BETAT*BETAT}",
                       abs(bT_f - BETAT) <= TOL_BETAT * BETAT))
        checks.append((run, "fit Tmd0", tmd_f, f"{TMD0}+-{TOL_TMD0}",
                       abs(tmd_f - TMD0) <= TOL_TMD0))
        checks.append((run, "fit exponent q", q_f, f"{Q}+-{TOL_Q}",
                       abs(q_f - Q) <= TOL_Q))
        fit_results[run] = dict(betaT=bT_f, Tmd0=tmd_f, q=q_f)
    else:
        m = slice(1, -1)
        c2, c1, c0 = np.polyfit(last["s"][m], b_rec[m], 2)
        # b(s) = -betaT*(1-Tmd0-slope*s)^2 + betaS*s expands to
        #   c0 = -betaT*(1-Tmd0)^2, c1 = betaS + 2*betaT*(1-Tmd0)*slope,
        #   c2 = -betaT*slope^2;  a slope-sign bug moves c1 by 1.2.
        c0_ex = -BETAT * (1.0 - TMD0) ** 2
        c1_ex = BETAS + 2.0 * BETAT * (1.0 - TMD0) * SLOPE
        c2_ex = -BETAT * SLOPE ** 2
        for cf, ce, name in [(c0, c0_ex, "c0"), (c1, c1_ex, "c1 (slope-sign)"),
                             (c2, c2_ex, "c2 (betaT*slope^2)")]:
            checks.append((run, f"fit {name}", cf, f"{ce}+-{TOL_C_RUNB}",
                           abs(cf - ce) <= TOL_C_RUNB))
        # implied EOS constants
        slope_mag = np.sqrt(max(-c2, 0.0) / BETAT)
        betaS_f = c1 - 2.0 * BETAT * (1.0 - TMD0) * np.sign(SLOPE) * slope_mag
        fit_results[run] = dict(c0=c0, c1=c1, c2=c2,
                                slope_mag=slope_mag, betaS=betaS_f)

    for i in range(len(x)):
        rows_profiles.append([run, x[i], last["th"][i], last["s"][i],
                              last["F"][i], w[i], w_bvp[i], b_in[i],
                              b_rec[i] if 0 < i < len(x) - 1 else np.nan])

# 7. optional MPI rank-count check
mpi_dir = "runA_theta_4rank"
if os.path.isdir(mpi_dir) and glob.glob(os.path.join(mpi_dir, "Data_*.h5")):
    f1 = sorted_outputs(mpi_dir)[-1]
    with h5py.File(f1, "r") as f:
        t1 = float(np.asarray(f["time"]).reshape(-1)[0])
    i4 = int(round(t1 / OUT_DT))
    f4 = os.path.join("runA_theta", f"Data_{i4}.h5")
    dmax = {}
    with h5py.File(f1, "r") as a, h5py.File(f4, "r") as b:
        for key in ["Conc/0", "Conc/1", "VOF/C_L", "u", "v", "w"]:
            dmax[key] = np.abs(a[key][:] - b[key][:]).max()
    worst = max(dmax.values())
    checks.append(("mpi", f"4-vs-1 rank max field diff at t={t1}",
                   worst, TOL_MPI, worst <= TOL_MPI))
    print("MPI check per field:", {k: f"{v:.3e}" for k, v in dmax.items()})
else:
    print("note: runA_theta_4rank/ not found - MPI check skipped")

# --- report -------------------------------------------------------------------
print(f"\n{'run':<12} {'criterion':<42} {'value':>13} {'threshold':>22} pass")
npass = 0
for run, crit, val, thr, ok in checks:
    print(f"{run:<12} {crit:<42} {val:>13.4e} {str(thr):>22} {'PASS' if ok else 'FAIL'}")
    npass += ok

with open("eos_results.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["run", "criterion", "value", "threshold", "pass"])
    for run, crit, val, thr, ok in checks:
        wcsv.writerow([run, crit, f"{val:.6e}", thr, int(ok)])
    for run in fit_results:
        for k, v in fit_results[run].items():
            wcsv.writerow([run, f"fitted_{k}", f"{v:.6f}", "", ""])

with open("eos_profiles.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["run", "x", "theta", "s", "F", "w", "w_bvp",
                   "b_input", "b_reconstructed"])
    wcsv.writerows(rows_profiles)

# --- figures ------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for col, run in enumerate(runs):
    first, prev, last = snap[run]
    x = last["x"]
    b_in = b_eos(last["th"], last["s"])
    w_bvp = solve_bvp(x, b_in)
    dx = x[1] - x[0]
    b_rec = np.full_like(last["w"], np.nan)
    b_rec[1:-1] = NU * (last["w"][2:] - 2 * last["w"][1:-1] + last["w"][:-2]) / dx**2

    ax = axes[0, col]
    ax.plot(x, last["w"], "b-", lw=1.5, label="w measured (steady)")
    ax.plot(x[::8], w_bvp[::8], "ko", ms=3, label="discrete BVP from inputs")
    ax.set_xlabel("x"); ax.set_ylabel("w")
    ax.set_title(f"{run}: steady profile (t = {last['t']:.2f})")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1, col]
    ax.plot(x, b_in, "r-", lw=1.5, label="b(theta,s) input EOS")
    ax.plot(x[1:-1:6], b_rec[1:-1:6], "ko", ms=3,
            label=r"$b_{rec}=\nu\,\delta_x^2 w$ measured")
    ax.set_xlabel("x"); ax.set_ylabel("b")
    ax.set_title("pointwise buoyancy reconstruction")
    ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("fig_eos_overview.png", dpi=150)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
first, prev, last = snap["runA_theta"]
x = last["x"]; dx = x[1] - x[0]
b_rec = NU * (last["w"][2:] - 2 * last["w"][1:-1] + last["w"][:-2]) / dx**2
th = np.linspace(0, 1, 200)
axes[0].plot(last["th"][1:-1], b_rec, "k.", ms=2.5, label="measured")
axes[0].plot(th, -BETAT * np.abs(th - TMD0) ** Q, "r-", lw=1,
             label=r"$-\beta_T|\theta-T_{md0}|^q$ input")
fr = fit_results["runA_theta"]
axes[0].axvline(TMD0, color="r", ls=":", lw=0.8)
axes[0].set_title(f"runA: b(theta); fit betaT={fr['betaT']:.4f}, "
                  f"Tmd0={fr['Tmd0']:.4f}, q={fr['q']:.4f}")
axes[0].set_xlabel(r"$\theta$"); axes[0].set_ylabel("b"); axes[0].legend()
axes[0].grid(alpha=0.3)

first, prev, last = snap["runB_salt"]
x = last["x"]; dx = x[1] - x[0]
b_rec = NU * (last["w"][2:] - 2 * last["w"][1:-1] + last["w"][:-2]) / dx**2
ss = np.linspace(0, 1, 200)
axes[1].plot(last["s"][1:-1], b_rec, "k.", ms=2.5, label="measured")
axes[1].plot(ss, b_eos(1.0, ss), "r-", lw=1, label="input EOS")
axes[1].plot(ss, -BETAT * (1 - TMD0 + SLOPE * ss) ** 2 + BETAS * ss, "b--",
             lw=1, label="WRONG slope sign (ref)")
fr = fit_results["runB_salt"]
axes[1].set_title(f"runB: b(s); fit c=({fr['c0']:.4f}, {fr['c1']:.4f}, "
                  f"{fr['c2']:.4f})")
axes[1].set_xlabel("s"); axes[1].set_ylabel("b"); axes[1].legend()
axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("fig_eos_identification.png", dpi=150)

# error plots
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for run in runs:
    first, prev, last = snap[run]
    x = last["x"]; dx = x[1] - x[0]
    b_in = b_eos(last["th"], last["s"])
    w_bvp = solve_bvp(x, b_in)
    b_rec = NU * (last["w"][2:] - 2 * last["w"][1:-1] + last["w"][:-2]) / dx**2
    axes[0].semilogy(x, np.abs(last["w"] - w_bvp) + 1e-20, label=run)
    axes[1].semilogy(x[1:-1], np.abs(b_rec - b_in[1:-1]) + 1e-20, label=run)
axes[0].set_title("|w - w_BVP|"); axes[1].set_title("|b_rec - b_input|")
for ax in axes:
    ax.set_xlabel("x"); ax.legend(); ax.grid(alpha=0.3)
axes[0].axhline(TOL_BVP, color="r", ls="--", lw=0.8)
axes[1].axhline(TOL_POINTWISE, color="r", ls="--", lw=0.8)
fig.tight_layout()
fig.savefig("fig_eos_errors.png", dpi=150)

nfail = len(checks) - npass
print(f"\n{'='*30}\nEOS GATE: {npass}/{len(checks)} criteria passed"
      f" -> {'PASS' if nfail == 0 else 'FAIL'}\n{'='*30}")
sys.exit(0 if nfail == 0 else 1)
