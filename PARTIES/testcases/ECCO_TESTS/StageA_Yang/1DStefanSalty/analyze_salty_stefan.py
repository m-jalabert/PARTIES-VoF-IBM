#!/usr/bin/env python3
"""
1D binary (salty) Stefan validation gate — roadmap section A.3b.

Fresh ice melting into warm salty water, both phases semi-infinite, u = 0.
Checks against the exact coupled similarity solution (two-phase Neumann +
salt-dilution boundary condition + liquidus depression):

    X(t) = x_int + 2*lambda*sqrt(alpha*t),      alpha = 1/Pe_T
    theta_i = T_melt - liquidus_slope*s_i
    s_i     = s_inf / (1 + sqrt(pi)*lam_s*exp(lam_s^2)*(1+erf(lam_s)))
    lam*sqrt(pi)*e^{lam^2}/stefan
            = (theta_inf-theta_i)/(1+erf(lam)) + (theta_ice-theta_i)/erfc(lam)

Gate criteria (roadmap A.3b — all must hold, exit 0 = PASS):
  1. front law:      lambda_fit within 3% of exact
  2. theta collapse: max error <= 0.05 in BOTH phases (outside the band)
  3. s collapse:     max error <= 0.05 (liquid side, outside the band)
  4. interface:      mean s(X) within 15% of s_i; mean theta(X) within 0.03
  5. salt budget:    |Int s - Int s(0)| / Int s(0) <= 1e-8 for all outputs
  6. salt in ice:    max s beyond the band (x > X + 25 dx) <= 1e-8
  7. enthalpy:       |Int(theta + F/St) - initial| relative <= 1e-8
  8. quiescence:     max |u|,|v|,|w| < 1e-10

Writes stefan_salty_results.csv (gate summary), stefan_salty_timeseries.csv,
stefan_salty_profiles.csv, and fig_*.png plots.

Run from the testcase folder after the simulation: python3 analyze_salty_stefan.py
"""

import csv
import glob
import re
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.optimize import brentq, curve_fit
from scipy.special import erf, erfc

# --- parameters (must match parties.inp) -------------------------------------
ST       = 0.5          # stefan = cp*dT/L
LAM_LIQ  = 0.5          # liquidus_slope
T_MELT   = 0.0
PE_T     = 100.0
PE_S     = 1000.0
THETA_W  = 1.0          # liquid far field
THETA_I  = 0.0          # ice far field (theta_ice)
S_INF    = 1.0          # liquid far-field salinity
X_INT    = 0.75         # initial interface (vof_slab_x0 * Lx)
ALPHA    = 1.0 / PE_T
DSALT    = 1.0 / PE_S
LE       = PE_S / PE_T

T_FIT_MIN = 0.5         # fit window start (band relaxation transient)

# --- exact coupled similarity constants ---------------------------------------
def s_interface(lam):
    ls = lam * np.sqrt(LE)
    return S_INF / (1.0 + np.sqrt(np.pi) * ls * np.exp(ls**2) * (1.0 + erf(ls)))

def stefan_residual(lam, with_salt=True, with_solid=True):
    si = s_interface(lam) if with_salt else 0.0
    ti = T_MELT - LAM_LIQ * si if with_salt else T_MELT
    rhs = (THETA_W - ti) / (1.0 + erf(lam))
    if with_solid:
        rhs += (THETA_I - ti) / erfc(lam)
    return lam * np.sqrt(np.pi) * np.exp(lam**2) / ST - rhs

lam_ex   = brentq(stefan_residual, 1e-6, 1.0)
s_i_ex   = s_interface(lam_ex)
theta_i_ex = T_MELT - LAM_LIQ * s_i_ex
lam_s_ex = lam_ex * np.sqrt(LE)
lam_nosolid  = brentq(lambda l: stefan_residual(l, True, False), 1e-6, 1.0)
lam_saltfree = brentq(lambda l: stefan_residual(l, False, True), 1e-6, 1.0)

print("exact coupled constants:")
print(f"  lambda   = {lam_ex:.6f}   (no-solid-conduction: {lam_nosolid:.6f}, "
      f"salt-free: {lam_saltfree:.6f})")
print(f"  s_i      = {s_i_ex:.6f}")
print(f"  theta_i  = {theta_i_ex:.6f}")

# --- exact similarity profiles ------------------------------------------------
def theta_exact(eta):
    liq = THETA_W + (theta_i_ex - THETA_W) * (1.0 + erf(eta)) / (1.0 + erf(lam_ex))
    sol = THETA_I + (theta_i_ex - THETA_I) * erfc(eta) / erfc(lam_ex)
    return np.where(eta <= lam_ex, liq, sol)

def s_exact(eta_s):
    liq = S_INF + (s_i_ex - S_INF) * (1.0 + erf(eta_s)) / (1.0 + erf(lam_s_ex))
    return np.where(eta_s <= lam_s_ex, liq, 0.0)

# --- load outputs --------------------------------------------------------------
files = sorted(glob.glob("Data_*.h5"),
               key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))
if len(files) < 5:
    sys.exit(f"only {len(files)} outputs found — run the case first")

with h5py.File(files[0], "r") as f:
    xc = f["grid/xc"][:]
dx = xc[1] - xc[0]
nx = len(xc) - 1                     # exclude the east boundary node
x = xc[:nx]

rows = []
th_prof, s_prof, F_prof = [], [], []
t_all, X_all, umax_all = [], [], []
salt_all, enth_all, sice_all = [], [], []
thX_all, sX_all = [], []

for fn in files:
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        th = f["Conc/0"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        s = f["Conc/1"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        F = f["VOF/C_L"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        umax = max(np.abs(f["u"][:]).max(), np.abs(f["v"][:]).max(),
                   np.abs(f["w"][:]).max())

    # front position: linear interpolation of F = 0.5 (F decreases with x)
    idx = np.where((F[:-1] >= 0.5) & (F[1:] < 0.5))[0]
    X = np.nan
    if len(idx):
        i = idx[-1]
        X = x[i] + dx * (F[i] - 0.5) / (F[i] - F[i + 1])

    thX = np.interp(X, x, th) if np.isfinite(X) else np.nan
    sX = np.interp(X, x, s) if np.isfinite(X) else np.nan

    salt = np.sum(s) * dx
    enth = np.sum(th + F / ST) * dx
    in_ice = x > (X + 25.0 * dx) if np.isfinite(X) else x > X_INT + 25.0 * dx
    s_ice = np.abs(s[in_ice]).max() if in_ice.any() else 0.0

    t_all.append(t); X_all.append(X); umax_all.append(umax)
    salt_all.append(salt); enth_all.append(enth); sice_all.append(s_ice)
    thX_all.append(thX); sX_all.append(sX)
    th_prof.append(th); s_prof.append(s); F_prof.append(F)

t_all = np.array(t_all); X_all = np.array(X_all)
salt_all = np.array(salt_all); enth_all = np.array(enth_all)
thX_all = np.array(thX_all); sX_all = np.array(sX_all)
umax_all = np.array(umax_all); sice_all = np.array(sice_all)

# --- 1. front law fit ----------------------------------------------------------
mask = (t_all >= T_FIT_MIN) & np.isfinite(X_all)

def front(t, lam_fit, t0):
    return X_INT + 2.0 * lam_fit * np.sqrt(np.clip(ALPHA * (t - t0), 0.0, None))

(lam_fit, t0), _ = curve_fit(front, t_all[mask], X_all[mask], p0=[lam_ex, 0.0])
rms = np.sqrt(np.mean((front(t_all[mask], lam_fit, t0) - X_all[mask])**2))
err_lam = abs(lam_fit - lam_ex) / lam_ex

# origin-free instantaneous lambda (diagnostic): d/dt (X-x_int)^2 = 4 lam^2 alpha
lam_inst = np.full_like(t_all, np.nan)
ok = np.isfinite(X_all)
if ok.sum() > 3:
    dXsq = np.gradient((X_all[ok] - X_INT)**2, t_all[ok])
    lam_inst[ok] = np.sqrt(np.clip(dXsq / (4.0 * ALPHA), 0.0, None))

print(f"\n1. front law: lambda_fit = {lam_fit:.5f} (t0 = {t0:+.4f}), "
      f"rel.err = {err_lam*100:.2f}%  [exact {lam_ex:.5f}], fit rms = {rms:.2e} "
      f"(dx = {dx:.2e})")

# --- 2./3. similarity collapse ---------------------------------------------------
band_excl = 12.0 * dx
th_errs, s_errs, sim_times = [], [], []
for frac in (0.4, 0.6, 0.8, 1.0):
    n = int(frac * (len(files) - 1))
    t = t_all[n]
    if t <= max(t0, 0.0) + 0.1 or not np.isfinite(X_all[n]):
        continue
    eta = (x - X_INT) / (2.0 * np.sqrt(ALPHA * (t - t0)))
    eta_s = (x - X_INT) / (2.0 * np.sqrt(DSALT * (t - t0)))
    outside = np.abs(x - X_all[n]) > band_excl
    th_err = np.max(np.abs(th_prof[n][outside] - theta_exact(eta[outside])))
    liq = outside & (x < X_all[n])
    s_err = np.max(np.abs(s_prof[n][liq] - s_exact(eta_s[liq])))
    th_errs.append(th_err); s_errs.append(s_err); sim_times.append(t)
    print(f"2./3. t = {t:5.2f}: max|theta-exact| = {th_err:.3e} , "
          f"max|s-exact| = {s_err:.3e}")

# --- 4. interface values ----------------------------------------------------------
sX_mean = np.nanmean(sX_all[mask])
thX_mean = np.nanmean(thX_all[mask])
err_sX = abs(sX_mean - s_i_ex) / s_i_ex
err_thX = abs(thX_mean - theta_i_ex)
print(f"4. interface: mean s(X) = {sX_mean:.4f} [exact {s_i_ex:.4f}, "
      f"rel.err {err_sX*100:.1f}%] ; mean theta(X) = {thX_mean:.4f} "
      f"[exact {theta_i_ex:.4f}, abs.err {err_thX:.4f}]")

# --- 5./7. conservation ------------------------------------------------------------
salt_drift = np.max(np.abs(salt_all - salt_all[0])) / abs(salt_all[0])
enth_drift = np.max(np.abs(enth_all - enth_all[0])) / abs(enth_all[0])
print(f"5. salt budget drift      : {salt_drift:.3e} relative")
print(f"7. enthalpy budget drift  : {enth_drift:.3e} relative")

# --- 6. salt in ice / 8. quiescence -------------------------------------------------
s_ice_max = sice_all.max()
umax = umax_all.max()
print(f"6. max salt in ice        : {s_ice_max:.3e}")
print(f"8. max |u|,|v|,|w|        : {umax:.3e}")

# --- verdict --------------------------------------------------------------------
checks = [
    ("1 front law lambda",      f"{lam_fit:.5f}", f"{lam_ex:.5f}",
     "rel.err <= 3%",   err_lam <= 0.03),
    ("2 theta similarity",      f"{max(th_errs):.3e}", "0",
     "max err <= 0.05", max(th_errs) <= 0.05),
    ("3 s similarity",          f"{max(s_errs):.3e}", "0",
     "max err <= 0.05", max(s_errs) <= 0.05),
    ("4a interface s(X)",       f"{sX_mean:.4f}", f"{s_i_ex:.4f}",
     "rel.err <= 15%",  err_sX <= 0.15),
    ("4b interface theta(X)",   f"{thX_mean:.4f}", f"{theta_i_ex:.4f}",
     "abs.err <= 0.03", err_thX <= 0.03),
    ("5 salt conservation",     f"{salt_drift:.3e}", "0",
     "<= 1e-8",         salt_drift <= 1e-8),
    ("6 salt in ice",           f"{s_ice_max:.3e}", "0",
     "<= 1e-8",         s_ice_max <= 1e-8),
    ("7 enthalpy conservation", f"{enth_drift:.3e}", "0",
     "<= 1e-8",         enth_drift <= 1e-8),
    ("8 quiescence",            f"{umax:.3e}", "0",
     "< 1e-10",         umax < 1e-10),
]

# --- CSV output -------------------------------------------------------------------
with open("stefan_salty_results.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["criterion", "measured", "reference", "tolerance", "pass"])
    for name, meas, ref, tol, ok_ in checks:
        wcsv.writerow([name, meas, ref, tol, "PASS" if ok_ else "FAIL"])
    wcsv.writerow([])
    wcsv.writerow(["fit lambda", f"{lam_fit:.6f}", f"{lam_ex:.6f}",
                   f"t0={t0:.4f}", f"rms={rms:.3e}"])
    wcsv.writerow(["lambda if solid conduction broken", f"{lam_nosolid:.6f}",
                   "", "", ""])
    wcsv.writerow(["lambda if salt coupling broken", f"{lam_saltfree:.6f}",
                   "", "", ""])

with open("stefan_salty_timeseries.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    wcsv.writerow(["t", "X_front", "lambda_inst", "theta_at_front", "s_at_front",
                   "umax", "salt_total", "enthalpy_total", "max_salt_in_ice"])
    for i in range(len(t_all)):
        wcsv.writerow([f"{t_all[i]:.6f}", f"{X_all[i]:.8f}", f"{lam_inst[i]:.6f}",
                       f"{thX_all[i]:.6f}", f"{sX_all[i]:.6f}",
                       f"{umax_all[i]:.3e}", f"{salt_all[i]:.12e}",
                       f"{enth_all[i]:.12e}", f"{sice_all[i]:.3e}"])

prof_idx = [int(f_ * (len(files) - 1)) for f_ in (0.0, 0.4, 0.7, 1.0)]
with open("stefan_salty_profiles.csv", "w", newline="") as fh:
    wcsv = csv.writer(fh)
    header = ["x"]
    for n in prof_idx:
        header += [f"theta_t{t_all[n]:.2f}", f"s_t{t_all[n]:.2f}",
                   f"F_t{t_all[n]:.2f}"]
    wcsv.writerow(header)
    for i in range(nx):
        row = [f"{x[i]:.8f}"]
        for n in prof_idx:
            row += [f"{th_prof[n][i]:.8e}", f"{s_prof[n][i]:.8e}",
                    f"{F_prof[n][i]:.8e}"]
        wcsv.writerow(row)

# --- plots -------------------------------------------------------------------------
tt = np.linspace(max(t0, 0.0) + 1e-6, t_all[-1], 300)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(t_all, X_all, "ko", ms=4, label="PARTIES (F = 0.5)")
ax.plot(tt, front(tt, lam_ex, t0), "r-", lw=1.5,
        label=fr"exact $\lambda$={lam_ex:.4f}")
ax.plot(tt, front(tt, lam_nosolid, t0), "b--", lw=1,
        label=fr"no solid conduction $\lambda$={lam_nosolid:.4f}")
ax.plot(tt, front(tt, lam_saltfree, t0), "g--", lw=1,
        label=fr"no salt coupling $\lambda$={lam_saltfree:.4f}")
ax.set_xlabel("t"); ax.set_ylabel("X(t)")
ax.set_title(fr"Front law: $\lambda_{{fit}}$ = {lam_fit:.4f} "
             fr"({err_lam*100:.2f}% err)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_front_law.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5))
eta_plot = np.linspace(-2.5, 2.5, 500)
ax.plot(eta_plot, theta_exact(eta_plot), "r-", lw=2, label="exact")
for frac in (0.4, 0.6, 0.8, 1.0):
    n = int(frac * (len(files) - 1))
    t = t_all[n]
    if t <= max(t0, 0.0) + 0.1:
        continue
    eta = (x - X_INT) / (2.0 * np.sqrt(ALPHA * (t - t0)))
    sel = np.abs(eta) < 2.5
    ax.plot(eta[sel], th_prof[n][sel], ".", ms=2, label=f"t = {t:.2f}")
ax.axvline(lam_ex, color="k", ls=":", lw=0.8)
ax.axhline(theta_i_ex, color="k", ls=":", lw=0.8)
ax.set_xlabel(r"$\eta = (x-x_{int})/2\sqrt{\alpha (t-t_0)}$")
ax.set_ylabel(r"$\theta$")
ax.set_title("Temperature similarity collapse (both phases)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_similarity_theta.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5))
etas_plot = np.linspace(-4.0, lam_s_ex, 500)
ax.plot(etas_plot, s_exact(etas_plot), "r-", lw=2, label="exact")
for frac in (0.4, 0.6, 0.8, 1.0):
    n = int(frac * (len(files) - 1))
    t = t_all[n]
    if t <= max(t0, 0.0) + 0.1:
        continue
    eta_s = (x - X_INT) / (2.0 * np.sqrt(DSALT * (t - t0)))
    sel = (eta_s > -4.0) & (x < X_all[n] + 3 * dx)
    ax.plot(eta_s[sel], s_prof[n][sel], ".", ms=2, label=f"t = {t:.2f}")
ax.axvline(lam_s_ex, color="k", ls=":", lw=0.8)
ax.axhline(s_i_ex, color="k", ls=":", lw=0.8)
ax.set_xlabel(r"$\eta_s = (x-x_{int})/2\sqrt{D (t-t_0)}$")
ax.set_ylabel("s")
ax.set_title("Salinity similarity collapse (meltwater dilution)")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("fig_similarity_salt.png", dpi=150); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].plot(t_all[mask], sX_all[mask], "bo-", ms=3, label="s(X)")
axes[0].axhline(s_i_ex, color="r", ls="--", label=f"exact s_i = {s_i_ex:.4f}")
axes[0].set_xlabel("t"); axes[0].set_ylabel("s at front"); axes[0].legend()
axes[0].grid(alpha=0.3)
axes[1].plot(t_all[mask], thX_all[mask], "bo-", ms=3, label=r"$\theta$(X)")
axes[1].axhline(theta_i_ex, color="r", ls="--",
                label=fr"exact $\theta_i$ = {theta_i_ex:.4f}")
axes[1].set_xlabel("t"); axes[1].set_ylabel(r"$\theta$ at front")
axes[1].legend(); axes[1].grid(alpha=0.3)
fig.suptitle("Interface values (liquidus + dilution pinning)")
fig.tight_layout(); fig.savefig("fig_interface_values.png", dpi=150); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].semilogy(t_all[1:],
                 np.abs(salt_all[1:] - salt_all[0]) / abs(salt_all[0]) + 1e-18,
                 "b.-", label="salt")
axes[0].semilogy(t_all[1:],
                 np.abs(enth_all[1:] - enth_all[0]) / abs(enth_all[0]) + 1e-18,
                 "r.-", label="enthalpy")
axes[0].axhline(1e-8, color="k", ls=":", label="tolerance")
axes[0].set_xlabel("t"); axes[0].set_ylabel("relative drift")
axes[0].legend(); axes[0].grid(alpha=0.3)
axes[0].set_title("Global conservation (all-no-flux box)")
axes[1].semilogy(t_all, umax_all + 1e-18, "k.-")
axes[1].axhline(1e-10, color="k", ls=":")
axes[1].set_xlabel("t"); axes[1].set_ylabel("max |u|,|v|,|w|")
axes[1].grid(alpha=0.3); axes[1].set_title("Quiescence")
fig.tight_layout(); fig.savefig("fig_conservation.png", dpi=150); plt.close(fig)

fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
for n in prof_idx:
    axes[0].plot(x, th_prof[n], label=f"t = {t_all[n]:.2f}")
    axes[1].plot(x, s_prof[n])
    axes[2].plot(x, F_prof[n])
axes[0].set_ylabel(r"$\theta$"); axes[0].legend(); axes[0].grid(alpha=0.3)
axes[1].set_ylabel("s"); axes[1].grid(alpha=0.3)
axes[2].set_ylabel("F"); axes[2].set_xlabel("x"); axes[2].grid(alpha=0.3)
axes[0].set_title("Profiles (y,z-averaged)")
fig.tight_layout(); fig.savefig("fig_overview.png", dpi=150); plt.close(fig)

# --- final verdict ------------------------------------------------------------------
print("\n--- gate summary ---")
all_ok = True
for name, meas, ref, tol, ok_ in checks:
    print(f"  {name:26s} {meas:>12s}  (ref {ref}, {tol})  "
          f"{'PASS' if ok_ else 'FAIL'}")
    all_ok &= ok_
print("\nPASS" if all_ok else "\nFAIL")
sys.exit(0 if all_ok else 1)
