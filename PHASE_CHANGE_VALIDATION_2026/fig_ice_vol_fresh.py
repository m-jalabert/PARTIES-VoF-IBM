#!/usr/bin/env python3
"""Ice volume V(t)/V0 vs time for the freshwater case (Sm=0) at every
available resolution, against Yang et al.'s own digitized Sm=0 curve.

N=216,288,360,512,720: computed from surviving field snapshots under
YangFresh_07272026/N*/Data_*.h5 -- the resolution ladder. Cn_test (a
Cahn number sensitivity check, not a resolution rung) and prov1440 (only
8 early frames, not a full curve) are both excluded.

N=1440: the real converged production run, computed the same way from
YangFresh_07272026's sibling directory Yang_Sm0_07182026/Data_*.h5 (101
full field snapshots, confirmed NXM=NYM=1440 in its parties.inp) -- not
the YangFresh_07272026/prov1440 stub. All six use the same vol_series
integral, so they're directly comparable.

Yang reference: a directly digitized curve, not a derived point --
ECCO_TESTS/Yang_production/yang_fig3a_reference.csv, Sm=0 rows, digitized
from the arXiv vector figure 3(a) [=JFM fig. 4(a)] at 900 dpi (see
YANG_REFERENCE_CORRECTION_REPORT.md). This is the strongest of the Yang
comparison targets in the whole campaign, unlike the derived Sm=5/dSv=5
point used in fig_ice_vol.py.

All sources stop at t=200, the published comparison window.
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


def th(t, v):
    for i in range(1, len(v)):
        if v[i - 1] > 0.5 >= v[i]:
            return t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1])
    return float("nan")


def fig_ice_vol_fresh():
    resolutions = [216, 288, 360, 512, 720]
    ladder = {n: vol_series(SCR / f"YangFresh_07272026/N{n}") for n in resolutions}
    t1440, v1440 = vol_series(SCR / "Yang_Sm0_07182026")

    raw = [ln for ln in open(TC / "Yang_production/yang_fig3a_reference.csv")
           if not ln.startswith("#")]
    yrows = [r for r in csv.DictReader(raw) if r["Sm_g_per_kg"] == "0"]
    ty = np.array([float(r["time_t_ff"]) for r in yrows])
    vy = np.array([float(r["ice_volume_ratio"]) for r in yrows])

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(ty, vy, ls=LS[1], color=GREY, lw=1.8,
            label=r"Yang et al., $S_m=0$")
    for k, n in enumerate(resolutions):
        t, v = ladder[n]
        ax.plot(t, v, ls=LS[k % len(LS)], color=C[k % len(C)],
                marker=MK[k % len(MK)], ms=3, markevery=max(1, len(t) // 15),
                label=rf"$N={n}^2$")
    ax.plot(t1440, v1440, ls=LS[len(resolutions) % len(LS)],
            color=C[len(resolutions) % len(C)],
            marker=MK[len(resolutions) % len(MK)], ms=3,
            markevery=max(1, len(t1440) // 15), label=r"$N=1440^2$")
    ax.axhline(0.5, color=GREY, ls=":", lw=0.9)
    ax.set_xlim(0, 200); ax.set_ylim(0.15, 1.02)
    ax.set_xlabel(r"$t$ (free-fall units)"); ax.set_ylabel(r"$V(t)/V_0$")
    ax.set_title(r"$S_m=0$ (freshwater) -- resolution dependence", loc="left")
    ax.legend(frameon=False, loc="upper right", fontsize=7)
    ax.set_box_aspect(1)
    style(ax)

    fig.tight_layout()
    fig.savefig(FIG / "ice_vol_fresh.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(FIG / "ice_vol_fresh.png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    parts = [f"t_half({n})={th(*ladder[n]):.2f}" for n in resolutions]
    parts.append(f"t_half(1440)={th(t1440, v1440):.2f}")
    parts.append(f"t_half(Yang,Sm0)={th(ty, vy):.2f}")
    print(f"  ice_vol_fresh.{{pdf,png}}   " + "  ".join(parts))


if __name__ == "__main__":
    fig_ice_vol_fresh()
