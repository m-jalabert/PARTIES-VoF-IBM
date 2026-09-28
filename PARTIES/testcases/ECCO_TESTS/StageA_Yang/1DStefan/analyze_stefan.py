#!/usr/bin/env python3
"""
1D Stefan (Neumann) validation gate — roadmap section A.3.

Checks, against the exact one-phase Neumann solution with
lambda*exp(lambda^2)*erf(lambda) = St/sqrt(pi):

  1. Front law      X(t) = 2*lambda*sqrt(alpha*(t - t0)), alpha = 1/Pe_T,
                    fitted with a virtual time origin t0 (the tanh IC is not
                    the similarity profile, so the front needs a fit origin).
  2. Similarity     theta(x,t) collapses on 1 - erf(eta)/erf(lambda),
                    eta = x/(2*sqrt(alpha*(t - t0))).
  3. Quiescence     max |u|,|v|,|w| < 1e-10 for all outputs.
  4. Enthalpy       d/dt Int(theta + F/St) dV = hot-wall conductive influx.

Run from the testcase folder after the simulation: python3 analyze_stefan.py
"""

import glob
import re
import sys

import h5py
import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.special import erf

# --- parameters (must match parties.inp) -----------------------------------
ST = 0.25            # stefan = cp*dT/L; St = stefan*(theta_w - theta_m)
PE = 10.0            # Pe_T
ALPHA = 1.0 / PE
THETA_W = 1.0

# --- exact Neumann constant -------------------------------------------------
lam = brentq(lambda l: l * np.exp(l**2) * erf(l) - ST / np.sqrt(np.pi), 1e-6, 2.0)

# --- load outputs ------------------------------------------------------------
files = sorted(glob.glob("Data_*.h5"),
               key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))
if len(files) < 5:
    sys.exit(f"only {len(files)} outputs found — run the case first")

t_all, X_all, umax_all, E_all, Q_all, prof_all = [], [], [], [], [], []

with h5py.File(files[0], "r") as f:
    xc = f["grid/xc"][:]          # cell centers, length NX (last = wall node)
dx = xc[1] - xc[0]
nx_real = len(xc) - 1             # exclude the east boundary node

for fn in files:
    with h5py.File(fn, "r") as f:
        t = float(f["time"][()])
        th = f["Conc/0"][:]       # (NZ, NY, NX)
        F = f["VOF/C_L"][:]
        u = np.abs(f["u"][:]).max()
        v = np.abs(f["v"][:]).max()
        w = np.abs(f["w"][:]).max()

    thx = th[1:-1, 1:-1, :nx_real].mean(axis=(0, 1))   # y,z-averaged profiles
    Fx = F[1:-1, 1:-1, :nx_real].mean(axis=(0, 1))

    # front position: linear interpolation of F = 0.5 (F decreases with x)
    idx = np.where((Fx[:-1] >= 0.5) & (Fx[1:] < 0.5))[0]
    X = np.nan
    if len(idx):
        i = idx[-1]
        X = xc[i] + dx * (Fx[i] - 0.5) / (Fx[i] - Fx[i + 1])

    # enthalpy per unit cross-section and hot-wall conductive influx;
    # the ghost-node Dirichlet gives flux = (1/Pe)*(theta_w - theta_0)*2/dx
    E = np.sum(thx + Fx / ST) * dx
    Q = (1.0 / PE) * (THETA_W - thx[0]) * 2.0 / dx

    t_all.append(t); X_all.append(X); umax_all.append(max(u, v, w))
    E_all.append(E); Q_all.append(Q); prof_all.append(thx)

t_all = np.array(t_all); X_all = np.array(X_all)
E_all = np.array(E_all); Q_all = np.array(Q_all)

# --- 1. front law fit --------------------------------------------------------
mask = t_all > 2.0                       # skip the IC transient
def front(t, lam_fit, t0):
    return 2.0 * lam_fit * np.sqrt(ALPHA * (t - t0))

(lam_fit, t0), _ = curve_fit(front, t_all[mask], X_all[mask], p0=[lam, -0.05])
rms = np.sqrt(np.mean((front(t_all[mask], lam_fit, t0) - X_all[mask])**2))
err_lam = abs(lam_fit - lam) / lam

print(f"exact  lambda            : {lam:.4f}")
print(f"fitted lambda (t0={t0:+.3f}): {lam_fit:.4f}   rel.err = {err_lam*100:.2f}%")
print(f"front-fit RMS            : {rms:.2e}  (dx = {dx:.2e})")

# --- 2. similarity collapse --------------------------------------------------
def theta_exact(eta):
    th = 1.0 - erf(eta) / erf(lam)
    return np.where(eta <= lam, th, 0.0)

print("\nsimilarity profile (liquid side, eta < 0.9*lambda):")
sim_errs = []
for frac in (0.4, 0.6, 0.8, 1.0):
    n = int(frac * (len(files) - 1))
    t = t_all[n]
    if t <= t0:
        continue
    eta = xc[:nx_real] / (2.0 * np.sqrt(ALPHA * (t - t0)))
    sel = eta < 0.9 * lam
    e = np.max(np.abs(prof_all[n][sel] - theta_exact(eta[sel])))
    sim_errs.append(e)
    print(f"  t = {t:6.2f}: max|theta - exact| = {e:.3e}")

# --- 3. quiescence -----------------------------------------------------------
umax = np.nanmax(umax_all)
print(f"\nmax |u|,|v|,|w| over run  : {umax:.3e}")

# --- 4. enthalpy budget ------------------------------------------------------
# The wall flux decays like 1/sqrt(t), so trapezoid quadrature of the sparse
# outputs aliases the early transient.  Fit Q(t) = a/sqrt(t - b) on t >= 1
# (sub-1e-3 fit residual once self-similar) and close the budget exactly on
# the late window, where quadrature error is negligible.
mQ = t_all >= 1.0
(qa, qb), _ = curve_fit(lambda t, a, b: a / np.sqrt(t - b),
                        t_all[mQ], Q_all[mQ], p0=[0.35, -0.1])
q_rms = np.sqrt(np.mean((qa / np.sqrt(t_all[mQ] - qb) - Q_all[mQ])**2)) \
        / Q_all[mQ].mean()
t1 = t_all[mQ][0]
Qint_late = 2.0 * qa * (np.sqrt(t_all[-1] - qb) - np.sqrt(t1 - qb))
dE_late = E_all[-1] - E_all[np.argmin(np.abs(t_all - t1))]
closure = abs(Qint_late - dE_late) / max(abs(dE_late), 1e-30)
print(f"\nwall-flux fit rel.rms     : {q_rms:.2e}")
print(f"enthalpy gain  (t in [{t1:.0f},{t_all[-1]:.0f}]): {dE_late:.6e}")
print(f"wall influx    (t in [{t1:.0f},{t_all[-1]:.0f}]): {Qint_late:.6e}")
print(f"relative closure error    : {closure:.3e}")

# --- verdict -----------------------------------------------------------------
ok = (err_lam < 0.03) and (umax < 1e-10) and (closure < 1e-2) \
     and all(e < 0.05 for e in sim_errs)
print("\nPASS" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
