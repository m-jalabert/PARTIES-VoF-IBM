#!/usr/bin/env python3
"""
1D binary (salty) Stefan validation — roadmap section A.3b.

Three separate article-quality figures, each numerical vs. the exact coupled
similarity solution (two-phase Neumann + salt-dilution boundary condition +
liquidus depression):
  fig_front_article       — interface position X(t)
  fig_salinity_article    — salinity profile s(x) at the final output time
  fig_temperature_article — temperature profile theta(x) at the final output time

Run from the testcase folder after the simulation:
    python3 plot_front_salt_temp_article.py
"""

import glob
import re
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

from scipy.optimize import brentq, curve_fit
from scipy.special import erf, erfc

# --- physical parameters (must match parties.inp) -----------------------------
ST      = 0.5          # stefan = cp*dT/L
LAM_LIQ = 0.5          # liquidus_slope
T_MELT  = 0.0
PE_T    = 100.0
PE_S    = 1000.0
THETA_W = 1.0          # liquid far field
THETA_I = 0.0          # ice far field (theta_ice)
S_INF   = 1.0          # liquid far-field salinity
X_INT   = 0.75         # initial interface (vof_slab_x0 * Lx)
ALPHA   = 1.0 / PE_T
DSALT   = 1.0 / PE_S
LE      = PE_S / PE_T

T_FIT_MIN = 0.5        # fit window start (band relaxation transient)

# --- plot controls (article style) --------------------------------------------
DPI = 600
SAVE_PDF = True
SAVE_PNG = True
AXIS_LABEL_FONTSIZE = 27.0
TICK_LABEL_FONTSIZE = 21.0
PLOT_BOX_LINEWIDTH = 1.2
MAJOR_TICK_LENGTH = 6.0
MAJOR_TICK_WIDTH = 1.2
MINOR_TICK_LENGTH = 3.5
MINOR_TICK_WIDTH = 1.0
LEGEND_FONTSIZE = 17.0
MARKERSIZE = 8
PROFILE_MARKER_STRIDE = 10   # subsample the 768-point profiles so dots stay legible

# Shape-figure controls.
SHAPE_ANALYTICAL_COLOR = "0.25"
SHAPE_SOLID_COLOR = "0.0"
SHAPE_LINEWIDTH_NUMERICAL = 1.1
SHAPE_LINEWIDTH_ANALYTICAL = 1.4
SHAPE_LINEWIDTH_SOLID = 1.1


def style_axes(ax):
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE,
                    length=MAJOR_TICK_LENGTH, width=MAJOR_TICK_WIDTH,
                    direction="in", top=True, right=True)
    ax.tick_params(axis="both", which="minor",
                    length=MINOR_TICK_LENGTH, width=MINOR_TICK_WIDTH,
                    direction="in", top=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(PLOT_BOX_LINEWIDTH)


def save_fig(fig, stem):
    if SAVE_PNG:
        fig.savefig(f"{stem}.png", dpi=DPI, bbox_inches="tight")
    if SAVE_PDF:
        fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


# --- exact coupled similarity constants ---------------------------------------
def s_interface(lam):
    ls = lam * np.sqrt(LE)
    return S_INF / (1.0 + np.sqrt(np.pi) * ls * np.exp(ls**2) * (1.0 + erf(ls)))


def stefan_residual(lam):
    si = s_interface(lam)
    ti = T_MELT - LAM_LIQ * si
    rhs = (THETA_W - ti) / (1.0 + erf(lam)) + (THETA_I - ti) / erfc(lam)
    return lam * np.sqrt(np.pi) * np.exp(lam**2) / ST - rhs


lam_ex   = brentq(stefan_residual, 1e-6, 1.0)
s_i_ex   = s_interface(lam_ex)
theta_i_ex = T_MELT - LAM_LIQ * s_i_ex
lam_s_ex = lam_ex * np.sqrt(LE)


def theta_exact(eta):
    liq = THETA_W + (theta_i_ex - THETA_W) * (1.0 + erf(eta)) / (1.0 + erf(lam_ex))
    sol = THETA_I + (theta_i_ex - THETA_I) * erfc(eta) / erfc(lam_ex)
    return np.where(eta <= lam_ex, liq, sol)


def s_exact(eta_s):
    liq = S_INF + (s_i_ex - S_INF) * (1.0 + erf(eta_s)) / (1.0 + erf(lam_s_ex))
    return np.where(eta_s <= lam_s_ex, liq, 0.0)


# --- load outputs ---------------------------------------------------------------
files = sorted(glob.glob("Data_*.h5"),
               key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))
if len(files) < 5:
    sys.exit(f"only {len(files)} outputs found — run the case first")

with h5py.File(files[0], "r") as f:
    xc = f["grid/xc"][:]
dx = xc[1] - xc[0]
nx = len(xc) - 1                     # exclude the east boundary node
x = xc[:nx]

t_all, X_all = [], []
t_final, s_final, theta_final = None, None, None
for k, fn in enumerate(files):
    with h5py.File(fn, "r") as f:
        t = float(np.asarray(f["time"]).reshape(-1)[0])
        th = f["Conc/0"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        s = f["Conc/1"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))
        F = f["VOF/C_L"][:][1:-1, 1:-1, :nx].mean(axis=(0, 1))

    # front position: linear interpolation of F = 0.5 (F decreases with x)
    idx = np.where((F[:-1] >= 0.5) & (F[1:] < 0.5))[0]
    X = np.nan
    if len(idx):
        i = idx[-1]
        X = x[i] + dx * (F[i] - 0.5) / (F[i] - F[i + 1])

    t_all.append(t); X_all.append(X)
    if k == len(files) - 1:
        t_final, s_final, theta_final = t, s, th

t_all = np.array(t_all); X_all = np.array(X_all)

# --- front law: fit the time-origin shift t0 (lambda fixed at the exact value) --
mask = (t_all >= T_FIT_MIN) & np.isfinite(X_all)


def front(t, lam, t0):
    return X_INT + 2.0 * lam * np.sqrt(np.clip(ALPHA * (t - t0), 0.0, None))


(lam_fit, t0), _ = curve_fit(front, t_all[mask], X_all[mask], p0=[lam_ex, 0.0])
err_lam = abs(lam_fit - lam_ex) / lam_ex

print(f"exact lambda = {lam_ex:.6f}, fitted lambda = {lam_fit:.6f} "
      f"(rel.err {err_lam*100:.2f}%), fitted t0 = {t0:+.4f}")
print(f"final output: t = {t_final:.4f}")

# --- figure 1: front position ----------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
tt = np.linspace(max(t0, 0.0) + 1e-6, t_all[-1], 300)
ax.plot(tt, front(tt, lam_ex, t0), "-", color=SHAPE_ANALYTICAL_COLOR,
        lw=SHAPE_LINEWIDTH_ANALYTICAL, label="Analytical")
ax.plot(t_all, X_all, "o", color=SHAPE_SOLID_COLOR,
        ms=MARKERSIZE, mew=SHAPE_LINEWIDTH_NUMERICAL, label="PARTIES")
ax.set_xlabel(r"$t$", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylabel(r"$X(t)$", fontsize=AXIS_LABEL_FONTSIZE)
ax.legend(fontsize=LEGEND_FONTSIZE, frameon=False)
style_axes(ax)
fig.tight_layout()
save_fig(fig, "fig_front_article")

# --- figure 2: salinity profile at final time -------------------------------------
xx = np.linspace(x[0], x[-1], 500)
eta_s_plot = (xx - X_INT) / (2.0 * np.sqrt(DSALT * (t_final - t0)))
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(xx, s_exact(eta_s_plot), "-", color=SHAPE_ANALYTICAL_COLOR,
        lw=SHAPE_LINEWIDTH_ANALYTICAL, label="Analytical")
ax.plot(x[::PROFILE_MARKER_STRIDE], s_final[::PROFILE_MARKER_STRIDE], "o",
        color=SHAPE_SOLID_COLOR, ms=MARKERSIZE, mew=SHAPE_LINEWIDTH_NUMERICAL,
        linestyle="none", label="PARTIES")
ax.set_xlabel(r"$x$", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylabel(r"$s$", fontsize=AXIS_LABEL_FONTSIZE)
ax.legend(fontsize=LEGEND_FONTSIZE, frameon=False)
style_axes(ax)
fig.tight_layout()
save_fig(fig, "fig_salinity_article")

# --- figure 3: temperature profile at final time ----------------------------------
eta_plot = (xx - X_INT) / (2.0 * np.sqrt(ALPHA * (t_final - t0)))
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(xx, theta_exact(eta_plot), "-", color=SHAPE_ANALYTICAL_COLOR,
        lw=SHAPE_LINEWIDTH_ANALYTICAL, label="Analytical")
ax.plot(x[::PROFILE_MARKER_STRIDE], theta_final[::PROFILE_MARKER_STRIDE], "o",
        color=SHAPE_SOLID_COLOR, ms=MARKERSIZE, mew=SHAPE_LINEWIDTH_NUMERICAL,
        linestyle="none", label="PARTIES")
ax.set_xlabel(r"$x$", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylabel(r"$\theta$", fontsize=AXIS_LABEL_FONTSIZE)
ax.legend(fontsize=LEGEND_FONTSIZE, frameon=False)
style_axes(ax)
fig.tight_layout()
save_fig(fig, "fig_temperature_article")

print("wrote fig_front_article, fig_salinity_article, fig_temperature_article "
      f"({'png' if SAVE_PNG else ''}{' + ' if SAVE_PNG and SAVE_PDF else ''}"
      f"{'pdf' if SAVE_PDF else ''})")
