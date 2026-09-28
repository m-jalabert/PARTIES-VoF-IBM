#!/usr/bin/env python3
"""Coarse-vs-fine grid comparison for the Yang salty case (Sm=5, dSv=5).

Fine (top) = N=432, coarse (bottom) = N=288. The true N=1440 salty
production run has no surviving field snapshots (Data_*.h5 purged from
scratch; only scalar summaries remain -- see YANG_VALIDATION_SUMMARY.md
Sec. 5), so N=432 stands in as "fine": it is the finer of the two
resolutions for which full field data still exists, and it is the same
run already used as "PARTIES" in yang_vs_reference.pdf.

Snapshot time matches yang_vs_reference.pdf: t=252 (nearest available
frame in both runs' 176-frame output cadence). Same colour code as that
figure (Yang's: temperature=magma over [0,1], salinity=Blues over
[0, Sm+dSv/2] -- matching Yang fig. 3's own colourbar labels "0 S
Sm+dSv/2", not the narrower [Sm-dSv/2, Sm+dSv/2] used here previously),
same solid-masked-white convention.

Conventions (snapshot loader, colour code, masking) copied verbatim from
make_figures.py::_snap and ::_yang_pair to keep this figure visually
consistent with the rest of the article's Yang panels.
"""
from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
SCR = Path("/anvil/scratch/x-mjalabert")
DATA_RE = re.compile(r"Data_(\d+)\.h5$")

T_TARGET = 252.0  # same instant as yang_vs_reference.pdf

# Yang et al. colour code, identified by matching their colourbars.
YCM_T, YCM_S, YS_LIM = "magma", "Blues", (0.0, 1.5)

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.labelsize": 10,
    "axes.titlesize": 10, "legend.fontsize": 8, "xtick.labelsize": 9,
    "ytick.labelsize": 9, "axes.linewidth": 0.7, "figure.dpi": 200,
})


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


def _yang_pair(ax_t, ax_s, S, arrows=True):
    X, Y = S["x"], S["y"]
    ext = [X[0], X[-1], Y[0], Y[-1]]
    solid = S["VOF/C_L"] < 0.5
    th = np.ma.array(S["Conc/0"], mask=solid)
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


def vv0(S):
    return float((1 - S["VOF/C_L"]).sum() * (S["x"][1] - S["x"][0])
                 * (S["y"][1] - S["y"][0]) / 0.1)


def main():
    fine = _snap(SCR / "YangSalty_07282026/N432", T_TARGET)
    coarse = _snap(SCR / "YangSalty_07282026/N288", T_TARGET)

    fig = plt.figure(figsize=(7.0, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.06, wspace=0.05)
    ft, fs_ = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ct, cs_ = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    _yang_pair(ft, fs_, fine)
    _yang_pair(ct, cs_, coarse)

    ft.set_title(r"temperature $\theta$", fontsize=9)
    fs_.set_title(r"salinity $s$", fontsize=9)
    ft.set_ylabel(rf"$N=432$ (fine)" "\n" rf"$V/V_0={vv0(fine):.2f}$", fontsize=9)
    ct.set_ylabel(rf"$N=288$ (coarse)" "\n" rf"$V/V_0={vv0(coarse):.2f}$", fontsize=9)

    fig.suptitle(rf"$S_m=5$, $\Delta S_v=5$, $t={fine['t']:.0f}$", fontsize=9, y=0.99)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"yang_coarse_vs_fine.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"  yang_coarse_vs_fine.pdf/png   "
          f"fine(N432) t={fine['t']:.1f} V/V0={vv0(fine):.3f}  |  "
          f"coarse(N288) t={coarse['t']:.1f} V/V0={vv0(coarse):.3f}")


if __name__ == "__main__":
    main()
