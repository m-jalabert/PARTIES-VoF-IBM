#!/usr/bin/env python3
"""Ice volume V(t)/V0 vs time for the salty case (Sm=5, dSv=5) at N=288,
432, and 1440, against Yang et al.'s own n_x=288 curve.

288 and 432: computed here from surviving field snapshots
(YangSalty_07282026/{N288,N432}/Data_*.h5), same integral as
make_figures.py::vol_series.

1440: no field snapshots survive for this run (see
YANG_VALIDATION_SUMMARY.md Sec. 5 -- the salty 1440^2 production binary's
trajectory is gone). Its V(t)/V0 is read directly from
Yang_production/final_timeseries.csv (column VV0), which is the one
surviving artefact of that run and needs no recomputation. That CSV only
reaches t=200 (the production run's time_max), so the 1440^2 curve is
shorter than the other two, which run to t=350.

Yang reference: digitized from the paper's own JFM Figure 1(a) -- their
grid-independence test, run at exactly this Sm=5, dSv=5 case, with base
resolutions n_x=72..288 overlaid. The n_x=288 curve (their production
choice) was extracted by pixel classification against its legend swatch
color from the CC-BY figure embedded in the Cambridge Core HTML full text
(source archived at
docs/papers/yang2023_arxiv_source/yang_JFM_fig1_grid_independence.png).
See ECCO_TESTS/Yang_production/yang_fig1a_nx288_reference.csv for the full
digitization provenance and its cross-check (t_half=216.76 vs their own
panel-b value ~216.4, 0.2% agreement).
"""
from __future__ import annotations

import csv
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
SCR = Path("/anvil/scratch/x-mjalabert")
TC = Path("/home/x-mjalabert/PARTIES/PARTIES/testcases/ECCO_TESTS")
DATA_RE = re.compile(r"Data_(\d+)\.h5$")


def style(ax):
    ax.grid(True, color=GREY, alpha=0.18, linewidth=0.4)
    ax.set_axisbelow(True)


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


def fig_ice_vol():
    t288, v288 = vol_series(SCR / "YangSalty_07282026/N288")
    t432, v432 = vol_series(SCR / "YangSalty_07282026/N432")

    rows = list(csv.DictReader(open(SCR / "Yang_production/final_timeseries.csv")))
    t1440 = np.array([float(r["t"]) for r in rows])
    v1440 = np.array([float(r["VV0"]) for r in rows])

    yrows = list(csv.DictReader(
        ln for ln in open(TC / "Yang_production/yang_fig1a_nx288_reference.csv")
        if not ln.startswith("#")))
    ty = np.array([float(r["time_t_ff"]) for r in yrows])
    vy = np.array([float(r["V_V0"]) for r in yrows])

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(ty, vy, ls=LS[1], color=GREY, lw=1.8,
            label=r"Yang et al., $n_x=288$")
    ax.plot(t288, v288, ls=LS[0], color=C[0], marker=MK[0], ms=3, markevery=8,
            label=r"$N=288^2$")
    ax.plot(t432, v432, ls=LS[2], color=C[1], marker=MK[1], ms=3, markevery=8,
            label=r"$N=432^2$")
    ax.plot(t1440, v1440, ls=LS[3], color=C[2], marker=MK[2], ms=3, markevery=8,
            label=r"$N=1440^2$")
    ax.axhline(0.5, color=GREY, ls=":", lw=0.9)
    ax.set_xlim(0, 300); ax.set_ylim(0.4, 1.02)
    ax.set_xlabel(r"$t$ (free-fall units)"); ax.set_ylabel(r"$V(t)/V_0$")
    ax.set_title(r"$S_m=5$, $\Delta S_v=5$ -- resolution dependence", loc="left")
    ax.legend(frameon=False, loc="upper right")
    ax.set_box_aspect(1)
    style(ax)

    fig.tight_layout()
    fig.savefig(FIG / "ice_vol.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(FIG / "ice_vol.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    def th(t, v):
        for i in range(1, len(v)):
            if v[i - 1] > 0.5 >= v[i]:
                return t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1])
        return float("nan")

    print(f"  ice_vol.{{pdf,png}}   "
          f"t_half(288)={th(t288, v288):.2f}  t_half(432)={th(t432, v432):.2f}  "
          f"t_half(1440)={th(t1440, v1440):.2f} (1440 run stops at t={t1440.max():.0f})  "
          f"t_half(Yang,288)={th(ty, vy):.2f}")


if __name__ == "__main__":
    fig_ice_vol()
