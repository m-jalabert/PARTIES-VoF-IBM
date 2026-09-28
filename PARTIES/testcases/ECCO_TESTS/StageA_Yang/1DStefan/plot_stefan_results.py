#!/usr/bin/env python3
"""
Presentation figures + CSV export for the 1D Stefan (Neumann) validation gate.

Reads the Data_*.h5 outputs in the current folder and produces:

  stefan_results.csv        time series: front position, instantaneous lambda,
                            enthalpy, wall flux, max velocity (+ exact values)
  stefan_profiles.csv       similarity profiles theta(eta) at several times
                            (+ the exact Neumann profile)
  fig_front_law.png         front X(t) vs the exact 2*lambda*sqrt(alpha*(t-t0))
  fig_lambda_inst.png       origin-free lambda(t) = sqrt(X*dX/dt/(2*alpha))
  fig_similarity.png        theta(eta) collapse onto 1 - erf(eta)/erf(lambda)
  fig_enthalpy_budget.png   enthalpy gain vs integrated hot-wall influx
  fig_overview.png          the four panels combined

Parameters must match parties.inp.  Run after the simulation:
    python3 plot_stefan_results.py
"""

import glob
import re

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq, curve_fit
from scipy.special import erf

# ---------------------------------------------------------------- parameters
ST = 0.25          # stefan = cp*dT/L  ->  St = stefan*(theta_w - theta_m)
PE = 10.0          # Pe_T
ALPHA = 1.0 / PE
THETA_W = 1.0
PROFILE_TIMES = (3.0, 6.0, 9.0, 12.0, 15.0)   # snapshots for the collapse plot

# exact Neumann constant:  lam*exp(lam^2)*erf(lam) = St/sqrt(pi)
LAM = brentq(lambda l: l * np.exp(l**2) * erf(l) - ST / np.sqrt(np.pi), 1e-6, 2.0)

# ---------------------------------------------------------------- load fields
files = sorted(glob.glob("Data_*.h5"),
               key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))

t_all, X_all, umax_all, E_all, Q_all, prof_all = [], [], [], [], [], []

with h5py.File(files[0], "r") as f:
    xc = f["grid/xc"][:]
dx = xc[1] - xc[0]
nx = len(xc) - 1                       # exclude the east boundary node

for fn in files:
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).ravel()[0])
        th = f["Conc/0"][:]
        F = f["VOF/C_L"][:]
        umax = max(np.abs(f[c][:]).max() for c in ("u", "v", "w"))

    thx = th[1:-1, 1:-1, :nx].mean(axis=(0, 1))
    Fx = F[1:-1, 1:-1, :nx].mean(axis=(0, 1))

    idx = np.where((Fx[:-1] >= 0.5) & (Fx[1:] < 0.5))[0]
    X = np.nan
    if len(idx):
        i = idx[-1]
        X = xc[i] + dx * (Fx[i] - 0.5) / (Fx[i] - Fx[i + 1])

    t_all.append(t)
    X_all.append(X)
    umax_all.append(umax)
    E_all.append(np.sum(thx + Fx / ST) * dx)                  # theta + F/St
    Q_all.append((1.0 / PE) * (THETA_W - thx[0]) * 2.0 / dx)  # ghost-Dirichlet flux
    prof_all.append(thx)

t_all = np.array(t_all)
X_all = np.array(X_all)
umax_all = np.array(umax_all)
E_all = np.array(E_all)
Q_all = np.array(Q_all)

# ------------------------------------------------------------- derived series
# front-law fit with virtual time origin (tanh IC != similarity profile)
fit_mask = t_all > 2.0
front = lambda t, lam, t0: 2.0 * lam * np.sqrt(ALPHA * (t - t0))
(lam_fit, t0), _ = curve_fit(front, t_all[fit_mask], X_all[fit_mask],
                             p0=[LAM, -0.05])

# origin-free instantaneous lambda: lam^2 = X * dX/dt / (2*alpha)
lam_inst = np.sqrt(np.clip(X_all * np.gradient(X_all, t_all) / (2 * ALPHA),
                           0.0, None))

# wall-flux fit a/sqrt(t-b) on t >= 1 and exact late-window influx integral
mQ = t_all >= 1.0
(qa, qb), _ = curve_fit(lambda t, a, b: a / np.sqrt(t - b),
                        t_all[mQ], Q_all[mQ], p0=[0.35, -0.1])
t1 = t_all[mQ][0]
Qint = np.where(t_all >= t1,
                2.0 * qa * (np.sqrt(np.maximum(t_all, t1) - qb)
                            - np.sqrt(t1 - qb)),
                np.nan)
dE = E_all - E_all[np.argmin(np.abs(t_all - t1))]
closure = abs(Qint[-1] - dE[-1]) / abs(dE[-1])

print(f"exact lambda  : {LAM:.4f}")
print(f"fitted lambda : {lam_fit:.4f}  (rel. err {abs(lam_fit-LAM)/LAM*100:.2f} %),"
      f"  t0 = {t0:+.3f}")
print(f"budget closure (t in [{t1:.0f},{t_all[-1]:.0f}]): {closure:.2e}")
print(f"max |u|,|v|,|w|: {umax_all.max():.2e}")

# ------------------------------------------------------------------ CSV no. 1
X_exact = front(t_all, LAM, t0)
header = ("time,front_X,front_X_exact,lambda_inst,lambda_exact,"
          "enthalpy,wall_flux,influx_int_late,dE_late,max_velocity")
table = np.column_stack([t_all, X_all, X_exact, lam_inst,
                         np.full_like(t_all, LAM), E_all, Q_all,
                         Qint, dE, umax_all])
np.savetxt("stefan_results.csv", table, delimiter=",", header=header,
           comments="", fmt="%.10e")

# ------------------------------------------------------------------ CSV no. 2
theta_exact = lambda eta: np.where(eta <= LAM,
                                   1.0 - erf(eta) / erf(LAM), 0.0)
eta_grid = np.linspace(0.0, 2.0 * LAM, 200)
cols, names = [eta_grid, theta_exact(eta_grid)], ["eta", "theta_exact"]
snap_idx = [int(np.argmin(np.abs(t_all - tt))) for tt in PROFILE_TIMES]
for n in snap_idx:
    eta = xc[:nx] / (2.0 * np.sqrt(ALPHA * (t_all[n] - t0)))
    cols.append(np.interp(eta_grid, eta, prof_all[n]))
    names.append(f"theta_t{t_all[n]:.0f}")
np.savetxt("stefan_profiles.csv", np.column_stack(cols), delimiter=",",
           header=",".join(names), comments="", fmt="%.10e")

# ------------------------------------------------------------------- plotting
plt.rcParams.update({
    "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13,
    "legend.fontsize": 11, "lines.linewidth": 1.8, "figure.dpi": 100,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_SIM, C_EX = "#1f77b4", "#d62728"


def panel_front(ax):
    tt = np.linspace(max(t0, 0.0) + 1e-6, t_all[-1], 400)
    ax.plot(tt, front(tt, LAM, t0), "-", color=C_EX,
            label=rf"Neumann: $2\lambda\sqrt{{\alpha(t-t_0)}}$, $\lambda={LAM:.4f}$")
    ax.plot(t_all, X_all, "o", ms=5, mfc="none", color=C_SIM,
            label=rf"PARTIES ($\lambda_{{\rm fit}}={lam_fit:.4f}$)")
    ax.set_xlabel(r"$t$"); ax.set_ylabel(r"front position $X(t)$")
    ax.set_title("Melt-front law"); ax.legend(frameon=False); ax.grid(alpha=0.3)


def panel_lambda(ax):
    m = t_all > 1.0
    ax.axhline(LAM, color=C_EX, label=rf"exact $\lambda={LAM:.4f}$")
    ax.plot(t_all[m], lam_inst[m], "o-", ms=4, color=C_SIM,
            label=r"$\sqrt{X\,\dot X/2\alpha}$ (origin-free)")
    ax.set_ylim(LAM * 0.96, LAM * 1.04)
    ax.set_xlabel(r"$t$"); ax.set_ylabel(r"$\lambda$")
    ax.set_title("Instantaneous front constant")
    ax.legend(frameon=False, loc="lower right"); ax.grid(alpha=0.3)


def panel_similarity(ax):
    ax.plot(eta_grid, theta_exact(eta_grid), "-", color=C_EX, zorder=5,
            label=r"$1-{\rm erf}(\eta)/{\rm erf}(\lambda)$")
    markers = "osd^v"
    for mk, n in zip(markers, snap_idx):
        eta = xc[:nx] / (2.0 * np.sqrt(ALPHA * (t_all[n] - t0)))
        sel = eta <= 2.0 * LAM
        ax.plot(eta[sel][::3], prof_all[n][sel][::3], mk, ms=4, mfc="none",
                label=rf"$t={t_all[n]:.0f}$")
    ax.axvline(LAM, color="gray", ls=":", lw=1)
    ax.text(LAM * 1.04, 0.38, r"$\eta=\lambda$ (front)", color="gray", fontsize=10)
    ax.set_xlabel(r"$\eta = x/2\sqrt{\alpha (t-t_0)}$"); ax.set_ylabel(r"$\theta$")
    ax.set_title("Self-similar temperature collapse")
    ax.legend(frameon=False); ax.grid(alpha=0.3)


def panel_budget(ax):
    m = t_all >= t1
    ax.plot(t_all[m], dE[m], "o", ms=5, mfc="none", color=C_SIM,
            label=r"enthalpy gain $\Delta\!\int(\theta + F/St)\,dV$")
    ax.plot(t_all[m], Qint[m], "-", color=C_EX,
            label="integrated hot-wall influx")
    ax.set_xlabel(r"$t$"); ax.set_ylabel(rf"budget since $t={t1:.0f}$")
    ax.set_title(f"Enthalpy conservation (closure {closure:.1e})")
    ax.legend(frameon=False); ax.grid(alpha=0.3)


for fname, fun in (("fig_front_law.png", panel_front),
                   ("fig_lambda_inst.png", panel_lambda),
                   ("fig_similarity.png", panel_similarity),
                   ("fig_enthalpy_budget.png", panel_budget)):
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    fun(ax)
    fig.savefig(fname)
    plt.close(fig)

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, fun in zip(axes.flat, (panel_front, panel_lambda,
                               panel_similarity, panel_budget)):
    fun(ax)
fig.suptitle("PARTIES 1D Stefan (Neumann) validation gate — "
             rf"$St={ST}$, $Pe_T={PE:.0f}$, $256\times4\times4$",
             fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("fig_overview.png")
plt.close(fig)

print("\nwritten: stefan_results.csv, stefan_profiles.csv, "
      "fig_front_law.png, fig_lambda_inst.png, fig_similarity.png, "
      "fig_enthalpy_budget.png, fig_overview.png")
