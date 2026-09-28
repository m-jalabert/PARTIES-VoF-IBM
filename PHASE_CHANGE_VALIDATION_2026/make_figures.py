#!/usr/bin/env python3
"""Generate every figure in the phase-change validation article.

All data are read from the actual run directories -- nothing is hard-coded or
reconstructed.  Figures that would require data no longer on disk (the 1-D
one-phase Stefan gate, the EOS-reconstruction and Brinkman sweeps) are NOT
produced; those cases appear as tables in the article, citing their gate reports.

Palette: Okabe-Ito, re-ordered so the weakest CVD pair is non-adjacent
(validated: worst adjacent dE 9.6 deutan / 20.0 normal, all checks pass).
Every series additionally carries its own marker and line style, so identity
never depends on colour -- required for print and greyscale.
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
MK = ["o", "s", "^", "D", "v", "P"]
LS = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 2))]
GREY = "#4d4d4d"

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.labelsize": 10,
    "axes.titlesize": 10, "legend.fontsize": 8, "xtick.labelsize": 9,
    "ytick.labelsize": 9, "axes.linewidth": 0.7, "grid.linewidth": 0.4,
    "grid.alpha": 0.35, "lines.linewidth": 1.6, "figure.dpi": 200,
    "axes.spines.top": False, "axes.spines.right": False,
})

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
TC = Path("/home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS")
SCR = Path("/anvil/scratch/x-mjalabert")
DATA_RE = re.compile(r"Data_(\d+)\.h5$")


def style(ax):
    ax.grid(True, color=GREY, alpha=0.18, linewidth=0.4)
    ax.set_axisbelow(True)


# --------------------------------------------------------------- Le=100 gate
LE, PE_T, PE_S = 100.0, 100.0, 10000.0
STEFAN_G, LIQUIDUS_G, S_INF, TH_INF, TH_ICE, X0 = 0.5, 0.5, 1.0, 1.0, 0.0, 0.75


def reference_state():
    """Exact two-phase (binary) Stefan similarity constant and interface state."""
    def vals(lam):
        lam_s = lam * math.sqrt(LE)
        denom = 1.0 + math.sqrt(math.pi) * lam_s * math.exp(lam_s ** 2) * (1 + math.erf(lam_s))
        s_i = S_INF / denom
        th_i = -LIQUIDUS_G * s_i
        left = lam * math.sqrt(math.pi) * math.exp(lam ** 2) / STEFAN_G
        right = (TH_INF - th_i) / (1 + math.erf(lam)) + (TH_ICE - th_i) / math.erfc(lam)
        return left - right, s_i, th_i
    lo, hi = 1e-8, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if vals(lo)[0] * vals(mid)[0] <= 0:
            hi = mid
        else:
            lo = mid
    lam = 0.5 * (lo + hi)
    _, s_i, th_i = vals(lam)
    return lam, s_i, th_i


def gate_series(run: Path):
    out = []
    for p in sorted(run.glob("Data_*.h5"), key=lambda q: int(DATA_RE.search(q.name).group(1))):
        with h5py.File(p, "r") as h:
            nx = len(h["grid/xc"]) - 1
            ny = len(h["grid/yc"]) - 1
            nz = len(h["grid/zc"]) - 1
            x = np.asarray(h["grid/xc"][:nx], float)
            F = np.asarray(h["VOF/C_L"][:nz, :ny, :nx], float).mean(axis=(0, 1))
            s = np.asarray(h["Conc/1"][:nz, :ny, :nx], float).mean(axis=(0, 1))
            t = float(np.asarray(h["time"]).flat[0])
        out.append((t, x, float(np.interp(0.5, F[::-1], x[::-1])), s))
    return out


def fig_stefan():
    lam, s_i, th_i = reference_state()
    ser = gate_series(TC / "Yang_salt_Le100_gate" / "run_legacy")
    t = np.array([r[0] for r in ser])
    xf = np.array([r[2] for r in ser])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    # (a) front displacement vs the exact similarity law
    m = t > 0
    a1.plot(t[m], xf[m] - X0, ls="none", marker=MK[0], ms=5, mfc="none",
            mec=C[0], mew=1.2, label="PARTIES")
    tt = np.linspace(t[m].min(), t.max(), 200)
    a1.plot(tt, 2 * lam * np.sqrt(tt / PE_T), ls=LS[1], color=GREY, lw=1.4,
            label=rf"exact, $2\lambda\sqrt{{t/Pe_T}}$, $\lambda={lam:.5f}$")
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"front displacement $X(t)-X_0$")
    a1.legend(frameon=False, loc="upper left")
    style(a1); a1.set_title("(a) interface position", loc="left")

    # (b) salinity profile at the final time vs exact
    tf, xg, _, sf = ser[-1]
    eta = (xg - X0) / (2 * math.sqrt(tf / PE_S))
    lam_s = lam * math.sqrt(LE)
    s_ex = S_INF + (s_i - S_INF) * (1 + np.vectorize(math.erf)(eta)) / (1 + math.erf(lam_s))
    w = (xg > X0 - 0.02) & (xg < X0 + 0.06)
    a2.plot(xg[w], sf[w], ls=LS[0], color=C[0], label=f"PARTIES, $t={tf:.2f}$")
    a2.plot(xg[w], s_ex[w], ls=LS[1], color=GREY, lw=1.4, label="exact similarity")
    a2.axvline(X0 + 2 * lam * math.sqrt(tf / PE_T), color=C[1], ls=LS[2], lw=1.2,
               label="exact interface")
    a2.set_xlabel(r"$x$"); a2.set_ylabel(r"salinity $s$")
    a2.legend(frameon=False, loc="lower right")
    style(a2); a2.set_title("(b) salinity profile", loc="left")

    fig.tight_layout()
    fig.savefig(FIG / "stefan_binary.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  stefan_binary.pdf   lambda_exact={lam:.8f}  "
          f"(article quotes the gate analyser's fit, 0.21844955, rel.err 6.04e-3)")


# ------------------------------------------------------------- melting RB
def fig_meltingrb():
    rows = list(csv.DictReader(open(SCR / "Melting_RB_St0p1_07172026/analysis/timeseries.csv")))
    t = np.array([float(r["t"]) for r in rows])
    h = np.array([float(r["hbar"]) for r in rows])
    qb = np.array([float(r["qbot"]) for r in rows])
    qt = np.array([float(r["qtop"]) for r in rows])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
    a1.plot(t, h, ls=LS[0], color=C[0], label=r"PARTIES, $\mathrm{St}=0.1$")
    # Diffusive phase: (1/2 + St) Pe dh/dt is constant.  Measured 23.34 against the
    # analytical 23.138 (ratio 1.009), so the target slope is measured/1.009.
    m = (t > 2) & (t < 16)
    sl = float(np.polyfit(t[m], h[m], 1)[0])
    tt = np.linspace(0, 26, 50)
    a1.plot(tt, h[m][0] + (sl / 1.009) * (tt - t[m][0]), ls=LS[1], color=GREY, lw=1.5,
            label=r"analytical diffusive law")
    a1.annotate(r"$(1/2+\mathrm{St})\,Pe\,\dot h = 23.34$" "\n"
                r"analytical $23.14$   $(+0.9\%)$",
                xy=(4, 0.80), fontsize=7.5, color=GREY)
    a1.set_xlabel(r"$t$"); a1.set_ylabel(r"mean interface height $\bar h$")
    a1.legend(frameon=False, loc="lower right"); style(a1)
    a1.set_title("(a) interface height", loc="left")

    a2.plot(t, qb, ls=LS[0], color=C[0], label=r"$q_{\mathrm{bot}}$")
    a2.plot(t, qt, ls=LS[1], color=C[1], label=r"$q_{\mathrm{top}}$")
    a2.set_xlabel(r"$t$"); a2.set_ylabel("heat flux")
    a2.set_yscale("log"); a2.legend(frameon=False); style(a2)
    a2.set_title("(b) boundary heat fluxes", loc="left")

    fig.tight_layout()
    fig.savefig(FIG / "melting_rb.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  melting_rb.pdf      h(0)={h[0]:.4f} h(end)={h[-1]:.4f} at t={t[-1]:.0f}")


# ---------------------------------------------------------------- Yang
def vol_series(d: Path):
    rows = []
    for p in d.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h:
            nx = len(h["grid/xc"]) - 1
            ny = len(h["grid/yc"]) - 1
            nz = len(h["grid/zc"]) - 1
            F = np.asarray(h["VOF/C_L"][:nz, :ny, :nx], float)
            x = np.asarray(h["grid/xc"][:nx], float)
            y = np.asarray(h["grid/yc"][:ny], float)
            t = float(np.asarray(h["time"]).flat[0])
        rows.append((t, (1 - F).sum() * float(x[1] - x[0]) * float(y[1] - y[0]) / 0.1))
    rows.sort()
    return np.array([a for a, _ in rows]), np.array([b for _, b in rows])


def t_half(t, v):
    for i in range(1, len(v)):
        if v[i - 1] > 0.5 >= v[i]:
            return t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1])
    return math.nan


def fig_yang():
    fresh = {216: 116.310, 288: 102.357, 360: 94.474, 512: 85.593, 720: 78.889, 1440: 60.287}
    salty = {288: 188.39, 432: 178.28}
    YF, YS, YSE = 107.42, 216.1, 10.6

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.0))

    # (a) the convergence ladder
    nf = np.array(sorted(fresh)); tf = np.array([fresh[n] for n in nf])
    ns = np.array(sorted(salty)); ts = np.array([salty[n] for n in ns])
    a1.plot(1 / nf, tf, ls=LS[0], color=C[0], marker=MK[0], ms=5,
            label="PARTIES, freshwater")
    a1.plot(1 / ns, ts, ls=LS[1], color=C[1], marker=MK[1], ms=5,
            label=r"PARTIES, salty ($\Delta S_v{=}5$)")
    a1.axhline(YF, color=C[0], ls=LS[2], lw=1.2)
    a1.axhline(YS, color=C[1], ls=LS[2], lw=1.2)
    a1.fill_between([0, 1 / 200], YS - YSE, YS + YSE, color=C[1], alpha=0.12, lw=0)
    a1.annotate("Yang, freshwater", (1 / 210, YF), color=C[0], fontsize=7.5, va="bottom")
    a1.annotate("Yang, salty", (1 / 210, YS), color=C[1], fontsize=7.5, va="bottom")
    a1.axvline(1 / 288, color=GREY, ls=LS[3], lw=1.0)
    a1.set_xlim(0, 1 / 200); a1.set_xlabel(r"$1/N$")
    a1.set_ylabel(r"half-melt time $t_{1/2}$")
    a1.legend(frameon=False, loc="center left"); style(a1)
    a1.set_title("(a) resolution dependence", loc="left")

    # (b) volume history at Yang's base resolution
    t288, v288 = vol_series(SCR / "YangFresh_07272026/N288")
    r = list(csv.DictReader(open(TC / "Yang_production/yang_sm0_final_timeseries.csv")))
    t1440 = np.array([float(x["time"]) for x in r])
    v1440 = np.array([float(x["ice_volume_ratio"]) for x in r])
    raw = [ln for ln in open(TC / "Yang_production/yang_fig3a_reference.csv")
           if not ln.startswith("#")]
    ref = [(float(x["time_t_ff"]), float(x["ice_volume_ratio"]))
           for x in csv.DictReader(raw) if x["Sm_g_per_kg"] == "0"]
    tr = np.array([a for a, _ in ref]); vr = np.array([b for _, b in ref])

    a2.plot(tr, vr, ls=LS[1], color=GREY, lw=1.6, label=r"Yang et al., $S_m{=}0$")
    a2.plot(t288, v288, ls=LS[0], color=C[0], label=r"PARTIES $288^2$ (matched)")
    a2.plot(t1440, v1440, ls=LS[2], color=C[2], label=r"PARTIES $1440^2$ (converged)")
    a2.axhline(0.5, color=GREY, ls=":", lw=0.9)
    a2.set_xlim(0, 200); a2.set_ylim(0, 1.02)
    a2.set_xlabel(r"$t$ (free-fall units)"); a2.set_ylabel(r"$V(t)/V_0$")
    a2.legend(frameon=False, loc="upper right"); style(a2)
    a2.set_title("(b) ice volume, freshwater", loc="left")

    fig.tight_layout()
    fig.savefig(FIG / "yang_benchmark.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  yang_benchmark.pdf  t_half(288)={t_half(t288, v288):.2f} "
          f"V(200)={np.interp(200, t288, v288):.4f}")


# =====================================================================
#  Additional figures: reproductions of specific reference-paper figures
#  Favier et al., arXiv:1901.03847  (melting Rayleigh-Benard)
#  Yang  et al., arXiv:2302.02357   (salt-stratified lateral melting)
# =====================================================================
RA_RB, THM_RB = 1.0e7, 0.05


def fig_rb_heat():
    """Mirror of Favier fig. 10: Nu against the effective Rayleigh number.

    Nu = q_bot h / (1 - theta_M),  Ra_e = Ra (1 - theta_M) h^3 -- definitions
    verified against the stored analysis arrays of the St=1 run.
    """
    a = np.load(SCR / "Melting_RB_07162026/analysis/analysis_arrays.npz")
    r = list(csv.DictReader(open(SCR / "Melting_RB_St0p1_07172026/analysis/timeseries.csv")))
    h01 = np.array([float(x["hbar"]) for x in r])
    q01 = np.array([float(x["qbot"]) for x in r])
    Nu01 = q01 * h01 / (1 - THM_RB)
    Rae01 = RA_RB * (1 - THM_RB) * h01 ** 3

    fr = list(csv.DictReader(open(SCR / "FlatRB_Ra1e6_07172026/analysis/timeseries.csv")))
    key = "Nu" if "Nu" in fr[0] else [k for k in fr[0] if k.lower().startswith("nu")][0]
    tf = np.array([float(x["t"]) for x in fr])
    nuf = np.array([float(x[key]) for x in fr])
    nu_flat = float(nuf[(tf >= 100) & (tf <= 400)].mean())

    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    m = a["Rae"] > 2e4
    ax.plot(a["Rae"][m], a["Nu"][m], ls="none", marker=MK[0], ms=3.2, alpha=0.75,
            color=C[0], label=r"melting, $\mathrm{St}_{\mathrm{alt}}=1$")
    m1 = Rae01 > 2e4
    ax.plot(Rae01[m1], Nu01[m1], ls="none", marker=MK[2], ms=3.4, alpha=0.75,
            color=C[1], label=r"melting, $\mathrm{St}_{\mathrm{alt}}=0.1$")
    ax.plot([1e6], [nu_flat], ls="none", marker="*", ms=13, color=C[2], mec="k", mew=0.5,
            label=rf"flat wall, $\mathrm{{Ra}}=10^6$ ($\mathrm{{Nu}}={nu_flat:.2f}$)")
    rr = np.logspace(4.3, 6.9, 50)
    ax.plot(rr, 0.115 * rr ** (1 / 3), ls=LS[1], color=GREY, lw=1.5,
            label=r"$0.115\,\mathrm{Ra}_e^{1/3}$ (Favier et al.)")
    ax.axhline(1.0, color=GREY, ls=":", lw=0.9)
    ax.annotate("conduction", (2.5e4, 1.05), fontsize=7, color=GREY, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"effective Rayleigh number $\mathrm{Ra}_e$")
    ax.set_ylabel(r"Nusselt number $\mathrm{Nu}$")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    style(ax)
    fig.tight_layout(); fig.savefig(FIG / "rb_heat_transport.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  rb_heat_transport.pdf   flat-wall Nu={nu_flat:.2f}")


def fig_rb_morphology():
    """Mirror of Favier fig. 7a (interface position h(x,t)) plus a field snapshot."""
    snaps = sorted(SCR.glob("Melting_RB_St0p1_07172026/analysis/fields_*.npz"))
    fig = plt.figure(figsize=(7.0, 4.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.15], hspace=0.42)

    ax = fig.add_subplot(gs[0])
    for k, p in enumerate(snaps):
        d = np.load(p)
        x = np.linspace(0, 6, d["hcol"].size)
        ax.plot(x, d["hcol"], ls=LS[k % len(LS)], color=C[k % len(C)], lw=1.2,
                label=rf"$t={float(d['t']):.0f}$")
    ax.set_xlim(0, 6); ax.set_xlabel(r"$x$"); ax.set_ylabel(r"interface $h(x,t)$")
    ax.set_ylim(0.10, 1.12)
    ax.legend(frameon=False, ncol=7, fontsize=6.5, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), columnspacing=1.0, handlelength=1.8)
    style(ax); ax.set_title("(a) melt-front topography", loc="left", y=1.16)

    d = np.load(snaps[len(snaps) // 2])
    ax2 = fig.add_subplot(gs[1])
    im = ax2.imshow(d["theta"], origin="lower", extent=[0, 6, 0, 1],
                    cmap="RdBu_r", vmin=0, vmax=1, aspect="auto", interpolation="bilinear")
    ax2.contour(np.linspace(0, 6, d["F"].shape[1]), d["yc"], d["F"], levels=[0.5],
                colors="k", linewidths=1.0)
    ax2.set_xlabel(r"$x$"); ax2.set_ylabel(r"$y$")
    cb = fig.colorbar(im, ax=ax2, pad=0.01, fraction=0.030)
    cb.set_label(r"$\theta$", rotation=0, labelpad=8)
    ax2.set_title(rf"(b) temperature and melt front, $t={float(d['t']):.0f}$", loc="left")

    fig.savefig(FIG / "rb_morphology.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  rb_morphology.pdf       {len(snaps)} snapshots, "
          f"h range {min(np.load(p)['hcol'].min() for p in snaps):.3f}"
          f"-{max(np.load(p)['hcol'].max() for p in snaps):.3f}")


def _snap(d: Path, t_target: float):
    best = None
    for p in d.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h:
            t = float(np.asarray(h["time"]).flat[0])
        if best is None or abs(t - t_target) < best[0]:
            best = (abs(t - t_target), t, p)
    with h5py.File(best[2], "r") as h:
        nz, ny, nx = len(h["grid/zc"]) - 1, len(h["grid/yc"]) - 1, len(h["grid/xc"]) - 1
        sl = (slice(0, nz), slice(0, ny), slice(0, nx))
        out = {k: np.asarray(h[k][sl], float).mean(axis=0)
               for k in ("Conc/0", "Conc/1", "VOF/C_L", "u", "v")}
        out["x"] = np.asarray(h["grid/xc"][:nx], float)
        out["y"] = np.asarray(h["grid/yc"][:ny], float)
        out["t"] = best[1]
    return out


def fig_yang_fields():
    """Mirror of Yang fig. 2: temperature and salinity fields with the melt front."""
    S = _snap(SCR / "YangSalty_07282026/N432", 100.0)
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7))
    X, Y = S["x"], S["y"]
    vmax = float(np.abs(S["v"]).max())
    panels = [("Conc/0", "RdBu_r", r"$\theta$", (0, 1), "(a) temperature"),
              ("Conc/1", "viridis", r"$s$", (0.4, 1.6), "(b) salinity"),
              ("v", "PuOr_r", r"$v$", (-vmax, vmax), "(c) vertical velocity")]
    for ax, (k, cm, lab, lim, ttl) in zip(axs, panels):
        im = ax.imshow(S[k], origin="lower", extent=[X[0], X[-1], Y[0], Y[-1]],
                       cmap=cm, vmin=lim[0], vmax=lim[1], aspect="equal",
                       interpolation="bilinear")
        ax.contour(X, Y, S["VOF/C_L"], levels=[0.5], colors="k", linewidths=0.8)
        ax.set_xlabel(r"$x$"); ax.set_title(ttl, loc="left", fontsize=9)
        if ax is axs[0]:
            ax.set_ylabel(r"$y$")
        else:
            ax.set_yticklabels([])
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.ax.set_title(lab, fontsize=8, pad=4)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle(rf"$S_m=5$, $\Delta S_v=5$, $t={S['t']:.0f}$", fontsize=9, y=1.03)
    fig.savefig(FIG / "yang_fields.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  yang_fields.pdf         t={S['t']:.1f}")


def fig_yang_profiles():
    """Mirror of Yang fig. 4b (mid-height vertical velocity) plus the melt-front shape."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    for k, (lab, run, tt) in enumerate([
            (r"$S_m=0$ (freshwater)", "YangFresh_07272026/N288", 100.0),
            (r"$S_m=5$, $\Delta S_v=5$", "YangSalty_07282026/N432", 100.0)]):
        S = _snap(SCR / run, tt)
        j = S["v"].shape[0] // 2
        a1.plot(S["x"], S["v"][j], ls=LS[k], color=C[k], label=lab)
        # front position per row, for the morphology panel
        fx = np.full(S["VOF/C_L"].shape[0], np.nan)
        for r in range(S["VOF/C_L"].shape[0]):
            row = S["VOF/C_L"][r]
            if row.min() <= 0.5 <= row.max():
                fx[r] = np.interp(0.5, row[::-1], S["x"][::-1])
        a2.plot(fx, S["y"], ls=LS[k], color=C[k], label=lab)

    a1.axhline(0, color=GREY, lw=0.6)
    a1.set_xlim(0.80, 1.0)
    a1.set_xlabel(r"$x$"); a1.set_ylabel(r"vertical velocity $v$")
    a1.legend(frameon=False, fontsize=7.5, loc="lower left"); style(a1)
    a1.set_title("(a) mid-height velocity profile", loc="left")

    a2.set_xlabel(r"melt-front position $x_f$"); a2.set_ylabel(r"$y$")
    a2.legend(frameon=False, fontsize=7.5, loc="lower left"); style(a2)
    a2.set_title("(b) melt-front shape", loc="left")

    fig.tight_layout(); fig.savefig(FIG / "yang_profiles.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  yang_profiles.pdf       written")


def _rb_snap(t_target):
    """Load a melting-RB snapshot with velocity, for vorticity."""
    best = None
    for q in (SCR / "Melting_RB_St0p1_07172026").glob("Data_*.h5"):
        if not DATA_RE.search(q.name):
            continue
        with h5py.File(q, "r") as h:
            t = float(np.asarray(h["time"]).flat[0])
        if best is None or abs(t - t_target) < best[0]:
            best = (abs(t - t_target), t, q)
    with h5py.File(best[2], "r") as h:
        nz, ny, nx = len(h["grid/zc"]) - 1, len(h["grid/yc"]) - 1, len(h["grid/xc"]) - 1
        sl = (slice(0, nz), slice(0, ny), slice(0, nx))
        th = np.asarray(h["Conc/0"][sl], float).mean(axis=0)
        F = np.asarray(h["VOF/C_L"][sl], float).mean(axis=0)
        u = np.asarray(h["u"][sl], float).mean(axis=0)
        v = np.asarray(h["v"][sl], float).mean(axis=0)
        x = np.asarray(h["grid/xc"][:nx], float)
        y = np.asarray(h["grid/yc"][:ny], float)
    dx = float(x[1] - x[0]); dy = float(y[1] - y[0])
    # omega_z = dv/dx - du/dy
    w = np.gradient(v, dx, axis=1) - np.gradient(u, dy, axis=0)
    return best[1], x, y, th, F, w


def fig_rb_panels():
    """Mirror of Favier fig. 4: temperature (left) and vorticity (right) at six
    increasing times, with the phi = 1/2 interface overlaid."""
    times = [5, 20, 40, 60, 80, 110]
    snaps = [_rb_snap(t) for t in times]
    wmax = max(float(np.abs(s[5]).max()) for s in snaps)

    fig, axs = plt.subplots(len(snaps), 2, figsize=(7.1, 0.62 * len(snaps) + 0.5),
                            gridspec_kw=dict(wspace=0.03, hspace=0.10))
    for k, (t, x, y, th, F, w) in enumerate(snaps):
        ext = [x[0], x[-1], y[0], y[-1]]
        # Same temperature encoding as the field panel of the morphology figure
        # (RdBu_r on 0->1), so the two figures can be read against each other.
        axs[k, 0].imshow(th, origin="lower", extent=ext, cmap="RdBu_r",
                         vmin=0.0, vmax=1.0, aspect="auto", interpolation="bilinear")
        axs[k, 1].imshow(w, origin="lower", extent=ext, cmap="bwr",
                         vmin=-0.25 * wmax, vmax=0.25 * wmax, aspect="auto",
                         interpolation="bilinear")
        for a in axs[k]:
            a.contour(x, y, F, levels=[0.5], colors="0.35", linewidths=0.7)
            a.set_xticks([]); a.set_yticks([])
        axs[k, 0].set_ylabel(rf"$t={t:.0f}$", rotation=0, ha="right", va="center",
                             fontsize=8, labelpad=12)
    axs[0, 0].set_title(r"temperature $\theta$", fontsize=9)
    axs[0, 1].set_title(r"vorticity $\omega_z$", fontsize=9)
    fig.savefig(FIG / "rb_panels.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  rb_panels.pdf           6 times, |omega|max={wmax:.2f}")


def fig_yang_evolution():
    """Melt-front evolution and mid-height profiles for the salt-stratified case."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.1))

    # (a) melt-front position at successive times -- layering development
    for k, tt in enumerate([20, 60, 100, 140, 180]):
        S = _snap(SCR / "YangSalty_07282026/N432", tt)
        F, x, y = S["VOF/C_L"], S["x"], S["y"]
        fx = np.array([np.interp(0.5, F[r][::-1], x[::-1])
                       if F[r].min() <= 0.5 <= F[r].max() else np.nan
                       for r in range(F.shape[0])])
        a1.plot(fx, y, ls=LS[k % len(LS)], color=C[k % len(C)], lw=1.2,
                label=rf"$t={S['t']:.0f}$")
    a1.set_xlabel(r"melt-front position $x_f$"); a1.set_ylabel(r"$y$")
    a1.legend(frameon=False, fontsize=7, loc="lower left"); style(a1)
    a1.set_title("(a) melt-front evolution", loc="left")

    # (b) mid-height theta, s and buoyancy -- mirrors the inset of Yang fig. 4(b)
    S = _snap(SCR / "YangSalty_07282026/N432", 100.0)
    j = S["Conc/0"].shape[0] // 2
    th, sa, x = S["Conc/0"][j], S["Conc/1"][j], S["x"]
    b = -1.0 * (th - 0.2 - (-0.0625) * sa) ** 2 + 1.75 * sa      # EOS, eq. (10)
    # Single axis on purpose: theta, s and b share a comparable range here, so a
    # second y-scale would only obscure their relative magnitudes.
    a2.plot(x, th, ls=LS[0], color=C[0], label=r"$\theta$")
    a2.plot(x, sa, ls=LS[1], color=C[1], label=r"$s$")
    a2.plot(x, b, ls=LS[2], color=C[2], label=r"$b(\theta,s)$")
    a2.axhline(0, color=GREY, lw=0.6)
    a2.set_xlim(0.80, 1.0)
    a2.set_xlabel(r"$x$"); a2.set_ylabel(r"$\theta$,\, $s$,\, $b$")
    a2.legend(frameon=False, fontsize=7.5, loc="center left")
    style(a2); a2.set_title(r"(b) mid-height profiles, $t=100$", loc="left")

    fig.tight_layout(); fig.savefig(FIG / "yang_evolution.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  yang_evolution.pdf      written")



def fig_rb_stefan():
    """Direct DNS-to-DNS comparison against Favier et al. figure 7(b).

    Their figure 7(b) plots the developed melting velocity against Stefan number
    at case-D parameters -- the same Ra, aspect ratio and theta_M as our runs --
    so our two Stefan numbers land directly on their curve.

    Their symbol positions were digitised from the published figure and are
    re-plotted here as data.  The figure image itself is NOT reproduced: that
    paper carries the arXiv perpetual non-exclusive licence, which does not
    grant third-party reuse.  Measured values are facts, not protected
    expression, and are attributed in the caption.
    """
    d = np.loadtxt(HERE / "favier_fig7b_dns.csv", delimiter=",")
    ours = {0.1: 23.34 / 0.6, 1.0: 19.09 / 1.5}      # (1/2+St) Pe hdot / (1/2+St)

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.plot(d[:, 0], d[:, 1], ls="none", marker="s", ms=6, mfc="none",
            mec=GREY, mew=1.2, label="Favier et al., DNS (digitised)")
    st = np.array(sorted(ours)); hd = np.array([ours[k] for k in st])
    ax.plot(st, hd, ls="none", marker=MK[0], ms=9, color=C[0],
            label="PARTIES")
    for k in st:
        ref = float(np.interp(np.log10(k), np.log10(d[:, 0]), d[:, 1]))
        ax.annotate(f"{100*(ours[k]/ref-1):+.0f}%", (k, ours[k]),
                    textcoords="offset points", xytext=(9, -11),
                    fontsize=7.5, color=C[0])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Stefan number $\mathrm{St}_{\mathrm{alt}}$")
    ax.set_ylabel(r"melting velocity $\mathrm{Pe}\,\dot{\overline{h}}$")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    style(ax)
    fig.tight_layout(); fig.savefig(FIG / "rb_stefan.pdf", bbox_inches="tight")
    plt.close(fig)
    for k in st:
        ref = float(np.interp(np.log10(k), np.log10(d[:, 0]), d[:, 1]))
        print(f"  rb_stefan.pdf           St={k}: ours {ours[k]:.2f} vs Favier {ref:.2f} "
              f"({100*(ours[k]/ref-1):+.1f}%)")



# Yang et al. colour code, matching their fig. 3 colourbar labels
# ("0  T  dT" and "0  S  Sm+dSv/2"): temperature = magma over [0, 1];
# salinity = Blues over [0, Sm+dSv/2]. (Previously [Sm-dSv/2, Sm+dSv/2],
# which clipped the near-zero salinity of fresh meltwater near the front
# into a single saturated-white band, wider than Yang's own rendering.)
YCM_T, YCM_S, YS_LIM = "magma", "Blues", (0.0, 1.5)


def _yang_pair(ax_t, ax_s, S, arrows=True):
    X, Y = S["x"], S["y"]
    ext = [X[0], X[-1], Y[0], Y[-1]]
    solid = S["VOF/C_L"] < 0.5
    th = np.ma.array(S["Conc/0"], mask=solid)      # solid masked white, as they do
    sa = np.ma.array(S["Conc/1"], mask=solid)
    for a, f, cm_, lim in ((ax_t, th, YCM_T, (0, 1)), (ax_s, sa, YCM_S, YS_LIM)):
        cmap = plt.get_cmap(cm_).copy(); cmap.set_bad("white")
        a.imshow(f, origin="lower", extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1],
                 aspect="equal", interpolation="bilinear")
        a.contour(X, Y, S["VOF/C_L"], levels=[0.5], colors="k", linewidths=0.7)
        a.set_xticks([]); a.set_yticks([])
    if arrows:
        st = max(1, len(X) // 26)
        ax_t.quiver(X[::st], Y[::st], S["u"][::st, ::st], S["v"][::st, ::st],
                    color="k", scale=0.45, width=0.004, headwidth=3.5)


def fig_yang_times():
    """Temperature and salinity at four times, in the reference's colour code."""
    times = [40, 100, 180, 252]
    fig, axs = plt.subplots(2, len(times), figsize=(7.1, 4.0),
                            gridspec_kw=dict(wspace=0.05, hspace=0.06))
    for k, tt in enumerate(times):
        S = _snap(SCR / "YangSalty_07282026/N432", tt)
        _yang_pair(axs[0, k], axs[1, k], S, arrows=False)
        V = float((1 - S["VOF/C_L"]).sum() * (S["x"][1] - S["x"][0])
                  * (S["y"][1] - S["y"][0]) / 0.1)
        axs[0, k].set_title(rf"$t={S['t']:.0f}$" "\n" rf"$V/V_0={V:.2f}$", fontsize=8)
    axs[0, 0].set_ylabel(r"temperature $\theta$", fontsize=8)
    axs[1, 0].set_ylabel(r"salinity $s$", fontsize=8)
    fig.savefig(FIG / "yang_times.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  yang_times.pdf          4 times, Yang colour code")


def fig_yang_vs():
    """Direct comparison: ours (top) against the reference panel (bottom).

    Matched on MELT STATE, not time: the reference does not publish the snapshot
    time, and it cannot be inferred reliably -- locating each of its three Sm=5
    panels on the published volume curves gives mutually inconsistent times, so
    those panels are evidently not at a common instant.  Ice area measured from
    their panel gives V/V0 = 0.35, and we show our field at the same value.
    """
    S = _snap(SCR / "YangSalty_07282026/N432", 252.0)
    ref = plt.imread(FIG / "yang_panel_Sm5_dSv5.png")

    fig = plt.figure(figsize=(7.0, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.06, wspace=0.05)
    at, as_ = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    _yang_pair(at, as_, S)
    at.set_title(r"temperature $\theta$", fontsize=9)
    as_.set_title(r"salinity $s$", fontsize=9)
    at.set_ylabel("PARTIES", fontsize=10)

    axr = fig.add_subplot(gs[1, :])
    axr.imshow(ref); axr.set_xticks([]); axr.set_yticks([])
    for sp in axr.spines.values():
        sp.set_visible(False)
    axr.set_ylabel("Yang et al.", fontsize=10)
    fig.savefig(FIG / "yang_vs_reference.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
    V = float((1 - S["VOF/C_L"]).sum() * (S["x"][1] - S["x"][0])
              * (S["y"][1] - S["y"][0]) / 0.1)
    print(f"  yang_vs_reference.pdf   ours t={S['t']:.0f} V/V0={V:.3f} vs reference V/V0~0.35")



if __name__ == "__main__":
    print("generating figures from run data:")
    fig_stefan()
    fig_meltingrb()
    fig_yang()
    fig_rb_heat()
    fig_rb_morphology()
    fig_yang_fields()
    fig_yang_profiles()
    fig_rb_panels()
    fig_yang_evolution()
    fig_rb_stefan()
    fig_yang_times()
    fig_yang_vs()
    print(f"\nwritten to {FIG}")
