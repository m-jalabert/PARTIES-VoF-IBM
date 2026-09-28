#!/usr/bin/env python3
"""
ICE_PENALIZATION validation gate - roadmap section A.3d.

Steady Brinkman-Poiseuille flow: water for x < 0.5 driven by a uniform
body force f = betaS = 2 (buoyancy with s = 1 in water, 0 in ice), ice slab
beyond, NO melting (stefan = 0).  The steady unidirectional flow w(x) obeys
the DISCRETE balance exactly (u = v = 0, dp = 0, advection = 0):

    (1/Re)*(w[i+1]-2w[i]+w[i-1])/dx^2 - (phi_i/tau)*w[i] + b_i = 0 ,
    phi = 1 - F  (ice fraction),  b = betaS * s ,

i.e. a Poiseuille profile in the liquid gap matched to an exponential
Brinkman layer of thickness delta = sqrt(nu*tau) inside the ice.  Sharp-
interface analytic reference (liquid gap X = 0.5, east wall at L = 1):

    liquid:  w = -(f/2nu) x^2 + A x ,  A = (fX/nu)(2 delta + X)/(2(delta+X))
    ice:     w = w_slip * exp(-(x-X)/delta) ,  w_slip ~ f X delta/(2 nu)

Three sub-runs sweep darcy_tau so that delta = 2, 4, 8 dx.

Gate criteria (exit 0 = PASS):
  1. steadiness:  max|w(t_end) - w(t_end-0.1)| <= 1e-5 * max|w|   (each tau)
  2. frozen state: F, s, theta drift <= 1e-6 (CH holds the slab; no melting)
  3. profile:     max|w - w_BVP| <= 1e-4 * max|w| with the BVP built from
                  the MEASURED F(x), s(x) (same discrete operator) (each tau)
  4. decay length: exponential fit deep in the ice (phi = 1) matches the
                  discrete Brinkman length to 2% and the continuum
                  delta = sqrt(nu*tau) to 5% (each tau) - a factor-2 error
                  in the Darcy diagonal would shift delta by 41%
  5. slip scaling: d ln(w_slip)/d ln(tau) = 0.5 +- 0.05 across the sweep
  6. ice rigidity: max|w| beyond 5 delta past the band <= 1% of max|w|
  7. transverse quiescence: max|u|, max|v| <= 1e-12 (each tau)
  8. [optional, if tau_4dx_4rank exists] MPI: 1-rank vs 4-rank fields agree
     to <= 1e-13

Writes penalization_results.csv, penalization_profiles.csv, fig_*.png.
Run from testcases/1DIcePenalization/: python3 analyze_penalization.py
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

# --- parameters (must match parties.inp) --------------------------------------
RE    = 10.0
NU    = 1.0 / RE
BETAS = 2.0          # force amplitude in the liquid
X_INT = 0.5          # nominal interface (vof_slab_x0 * Lx)
L     = 1.0
CN    = 0.0029296875
OUT_DT = 0.2         # output_time_interval (for the MPI-check index)

RUNS = {"tau_2dx": 6.103515625e-04,
        "tau_4dx": 2.44140625e-03,
        "tau_8dx": 9.765625e-03}

TOL_STEADY = 1.0e-5
TOL_FROZEN = 1.0e-6
TOL_BVP    = 1.0e-4      # relative on w
TOL_DDISC  = 0.02        # decay length vs discrete Brinkman length
TOL_DCONT  = 0.05        # decay length vs continuum sqrt(nu*tau)
TOL_SLOPE  = 0.05        # slip ~ tau^0.5
TOL_LEAK   = 0.01        # relative, beyond 5 delta past the band
TOL_UV     = 1.0e-12
TOL_MPI    = 1.0e-13


def load_snapshot(fn):
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        xc = f["grid/xc"][:]
        nx = len(xc) - 1
        th = f["Conc/0"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        s  = f["Conc/1"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        F  = f["VOF/C_L"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        w  = f["w"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        umax = np.abs(f["u"][:]).max()
        vmax = np.abs(f["v"][:]).max()
    return dict(t=t, x=xc[:nx], th=th, s=s, F=F, w=w, umax=umax, vmax=vmax)


def sorted_outputs(d):
    files = sorted(glob.glob(os.path.join(d, "Data_*.h5")),
                   key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))
    if len(files) < 3:
        sys.exit(f"{d}: only {len(files)} outputs - run the case first")
    return files


def solve_bvp(x, b, phi, tau):
    """nu*w'' - (phi/tau) w + b = 0, w = 0 at both walls (half-cell ghost)."""
    n = len(x)
    dx = x[1] - x[0]
    ab = np.zeros((3, n))
    ab[0, 1:] = NU / dx**2
    ab[2, :-1] = NU / dx**2
    ab[1, :] = -2.0 * NU / dx**2 - phi / tau
    ab[1, 0] -= NU / dx**2      # ghost = -w0
    ab[1, -1] -= NU / dx**2
    return solve_banded((1, 1), ab, -b)


def sharp_analytic(tau):
    """Sharp-interface reference: A (liquid slope) and slip w(X)."""
    delta = np.sqrt(NU * tau)
    f = BETAS
    X = X_INT
    A = (f * X / NU) * (2.0 * delta + X) / (2.0 * (delta + X))
    w_slip = X * (A - f * X / (2.0 * NU))
    return delta, A, w_slip


checks = []
rows_profiles = []
slip_meas, slip_bvp, taus = [], [], []
prof = {}

for run, tau in RUNS.items():
    files = sorted_outputs(run)
    first = load_snapshot(files[0])
    prev = load_snapshot(files[-2])
    last = load_snapshot(files[-1])
    prof[run] = last

    x, dx = last["x"], last["x"][1] - last["x"][0]
    wmax = np.abs(last["w"]).max()
    delta = np.sqrt(NU * tau)

    # 1. steadiness
    dsteady = np.abs(last["w"] - prev["w"]).max() / wmax
    checks.append((run, "steadiness d|w|/max|w|", dsteady, TOL_STEADY,
                   dsteady <= TOL_STEADY))

    # 2. frozen background state (no melting, frozen scalars, CH equilibrium)
    dF = np.abs(last["F"] - first["F"]).max()
    ds = np.abs(last["s"] - first["s"]).max()
    dth = np.abs(last["th"] - first["th"]).max()
    checks.append((run, "F drift (slab stationary)", dF, TOL_FROZEN, dF <= TOL_FROZEN))
    checks.append((run, "s drift", ds, TOL_FROZEN, ds <= TOL_FROZEN))
    checks.append((run, "theta drift", dth, TOL_FROZEN, dth <= TOL_FROZEN))

    # 3. w(x) against the BVP with the measured phi(x), b(x)
    phi = np.clip(1.0 - last["F"], 0.0, 1.0)
    b = BETAS * last["s"]
    w_bvp = solve_bvp(x, b, phi, tau)
    dbvp = np.abs(last["w"] - w_bvp).max() / wmax
    checks.append((run, "max|w - w_BVP|/max|w|", dbvp, TOL_BVP, dbvp <= TOL_BVP))

    # 4. exponential decay length deep in the ice (phi = 1 exactly)
    band = 6.0 * 2.0 * np.sqrt(2.0) * CN            # tanh band half-extent
    wabs = np.abs(last["w"])
    sel = (x > X_INT + band) & (phi > 1.0 - 1e-12) \
          & (wabs > 1e-11) & (wabs < 0.5 * np.abs(np.interp(X_INT, x, last["w"])))
    dfit = np.nan
    if sel.sum() >= 4:
        p = np.polyfit(x[sel], np.log(wabs[sel]), 1)
        dfit = -1.0 / p[0]
    delta_disc = dx / np.arccosh(1.0 + dx**2 / (2.0 * delta**2))
    e_disc = abs(dfit - delta_disc) / delta_disc
    e_cont = abs(dfit - delta) / delta
    checks.append((run, f"decay len vs discrete ({delta_disc:.5f})", e_disc,
                   TOL_DDISC, e_disc <= TOL_DDISC))
    checks.append((run, f"decay len vs sqrt(nu*tau) ({delta:.5f})", e_cont,
                   TOL_DCONT, e_cont <= TOL_DCONT))

    # 6. ice rigidity: leakage beyond 5 delta past the band
    far = x > X_INT + band + 5.0 * delta
    leak = wabs[far].max() / wmax if far.any() else 0.0
    checks.append((run, "ice leakage max|w|/max|w| beyond 5 delta", leak,
                   TOL_LEAK, leak <= TOL_LEAK))

    # 7. transverse quiescence
    checks.append((run, "max|u|", last["umax"], TOL_UV, last["umax"] <= TOL_UV))
    checks.append((run, "max|v|", last["vmax"], TOL_UV, last["vmax"] <= TOL_UV))

    # bookkeeping for the scaling law
    w_at_X = np.abs(np.interp(X_INT, x, last["w"]))
    slip_meas.append(w_at_X)
    slip_bvp.append(np.abs(np.interp(X_INT, x, w_bvp)))
    taus.append(tau)

    d_sharp, A_sharp, slip_sharp = sharp_analytic(tau)
    print(f"{run}: delta={delta:.5f} ({delta/dx:.1f} dx)  "
          f"w_max={wmax:.4f}  slip(X)={w_at_X:.5f} "
          f"(sharp {abs(slip_sharp):.5f}, bvp {slip_bvp[-1]:.5f})  "
          f"decay fit={dfit:.5f}  leak={leak:.2e}")

    for i in range(len(x)):
        rows_profiles.append([run, tau, x[i], last["F"][i], last["s"][i],
                              last["w"][i], w_bvp[i]])

# 5. slip ~ sqrt(tau) scaling across the sweep
lt, ls = np.log(np.array(taus)), np.log(np.array(slip_meas))
slope_fit = np.polyfit(lt, ls, 1)[0]
checks.append(("sweep", "d ln(slip)/d ln(tau)", slope_fit, f"0.5+-{TOL_SLOPE}",
               abs(slope_fit - 0.5) <= TOL_SLOPE))

# 8. optional MPI rank-count check
mpi_dir = "tau_4dx_4rank"
if os.path.isdir(mpi_dir) and glob.glob(os.path.join(mpi_dir, "Data_*.h5")):
    f1 = sorted_outputs(mpi_dir)[-1]
    with h5py.File(f1, "r") as f:
        t1 = float(np.asarray(f["time"]).reshape(-1)[0])
    i4 = int(round(t1 / OUT_DT))
    f4 = os.path.join("tau_4dx", f"Data_{i4}.h5")
    dmax = {}
    with h5py.File(f1, "r") as a, h5py.File(f4, "r") as b:
        for key in ["Conc/0", "Conc/1", "VOF/C_L", "u", "v", "w"]:
            dmax[key] = np.abs(a[key][:] - b[key][:]).max()
    worst = max(dmax.values())
    checks.append(("mpi", f"4-vs-1 rank max field diff at t={t1}", worst,
                   TOL_MPI, worst <= TOL_MPI))
    print("MPI check per field:", {k: f"{v:.3e}" for k, v in dmax.items()})
else:
    print("note: tau_4dx_4rank/ not found - MPI check skipped")

# --- report -------------------------------------------------------------------
print(f"\n{'run':<9} {'criterion':<46} {'value':>13} {'threshold':>13} pass")
npass = 0
for run, crit, val, thr, ok in checks:
    print(f"{run:<9} {crit:<46} {val:>13.4e} {str(thr):>13} {'PASS' if ok else 'FAIL'}")
    npass += ok

with open("penalization_results.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["run", "criterion", "value", "threshold", "pass"])
    for run, crit, val, thr, ok in checks:
        wcsv.writerow([run, crit, f"{val:.6e}", thr, int(ok)])
    wcsv.writerow([])
    wcsv.writerow(["run", "darcy_tau", "delta_theory", "slip_measured",
                   "slip_bvp", "slip_sharp"])
    for (run, tau), sm, sb in zip(RUNS.items(), slip_meas, slip_bvp):
        wcsv.writerow([run, tau, np.sqrt(NU * tau), f"{sm:.6e}",
                       f"{sb:.6e}", f"{abs(sharp_analytic(tau)[2]):.6e}"])
    wcsv.writerow(["sweep", "slip_scaling_exponent", f"{slope_fit:.4f}",
                   "0.5", ""])

with open("penalization_profiles.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["run", "darcy_tau", "x", "F", "s", "w", "w_bvp"])
    wcsv.writerows(rows_profiles)

# --- figures ------------------------------------------------------------------
colors = {"tau_2dx": "tab:blue", "tau_4dx": "tab:green", "tau_8dx": "tab:red"}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for run, tau in RUNS.items():
    last = prof[run]
    x = last["x"]
    phi = np.clip(1.0 - last["F"], 0.0, 1.0)
    w_bvp = solve_bvp(x, BETAS * last["s"], phi, tau)
    dlab = np.sqrt(NU * tau) / (x[1] - x[0])
    axes[0].plot(x, last["w"], color=colors[run], lw=1.5,
                 label=f"{run} (delta={dlab:.0f}dx)")
    axes[0].plot(x[::6], w_bvp[::6], "k.", ms=3)
    axes[1].semilogy(x, np.abs(last["w"]) + 1e-30, color=colors[run], lw=1.5)
    axes[1].semilogy(x[::4], np.abs(w_bvp[::4]) + 1e-30, "k.", ms=3)
axes[0].axvline(X_INT, color="gray", ls=":")
axes[0].set_xlabel("x"); axes[0].set_ylabel("w")
axes[0].set_title("steady Brinkman-Poiseuille profiles (dots: discrete BVP)")
axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].axvline(X_INT, color="gray", ls=":")
axes[1].set_xlabel("x"); axes[1].set_ylabel("|w|")
axes[1].set_ylim(1e-14, 1.0)
axes[1].set_title("exponential Brinkman decay in the ice")
axes[1].grid(alpha=0.3)
fig.tight_layout()
fig.savefig("fig_pen_profiles.png", dpi=150)

fig, ax = plt.subplots(figsize=(5.5, 4.5))
tt = np.array(taus)
ax.loglog(tt, slip_meas, "ko", label="measured w(X)")
ax.loglog(tt, [abs(sharp_analytic(t)[2]) for t in tt], "r-",
          label=r"sharp analytic $fX\delta/2\nu$ ($\propto\sqrt{\tau}$)")
ax.loglog(tt, slip_bvp, "b--", lw=1, label="diffuse-band BVP")
ax.set_xlabel(r"darcy_tau $\tau$"); ax.set_ylabel("interface slip |w(X)|")
ax.set_title(f"slip scaling: fitted exponent {slope_fit:.3f} (exact 1/2)")
ax.legend(); ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig("fig_pen_scaling.png", dpi=150)

nfail = len(checks) - npass
print(f"\n{'='*36}\nICE_PENALIZATION GATE: {npass}/{len(checks)} criteria "
      f"passed -> {'PASS' if nfail == 0 else 'FAIL'}\n{'='*36}")
sys.exit(0 if nfail == 0 else 1)
