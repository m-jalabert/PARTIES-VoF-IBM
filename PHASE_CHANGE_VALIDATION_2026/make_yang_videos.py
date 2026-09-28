#!/usr/bin/env python3
"""Temperature/salinity field videos for the Yang salty case (Sm=5, dSv=5),
one video per resolution: temperature (left) and salinity (right), side by
side, animated over all output times. No title, no legend, no colorbar --
just the fields, in Yang et al.'s own colour code, matching their fig. 3
colourbar labels ("0  T  dT", "0  S  Sm+dSv/2"): magma over [0,1] for
temperature, Blues over [0,1.5] (= [0, Sm+dSv/2] nondimensionalized by
Sm) for salinity. Solid (ice) is masked white; the melt front (C_L=0.5)
is drawn as a thin black contour.

Only N=288 and N=432 are made: the salty N=1440 production run has no
surviving field snapshots (Data_*.h5) to animate -- see
YANG_VALIDATION_SUMMARY.md Sec. 5. Scalar summaries survive for that run
but there is no field to show.
"""
from __future__ import annotations

import re
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
# The cluster's module-provided ffmpeg has no software H.264 encoder (only
# an unusable V4L2 hardware wrapper), which produced mpeg4-in-mp4 files
# that most standard players (QuickTime, Windows, browsers) can't open.
# imageio-ffmpeg's bundled static binary does have libx264, so use that.
import imageio_ffmpeg
matplotlib.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

HERE = Path(__file__).resolve().parent
VID = HERE / "videos"
VID.mkdir(exist_ok=True)
SCR = Path("/anvil/scratch/x-mjalabert")
DATA_RE = re.compile(r"Data_(\d+)\.h5$")

# Yang et al. colour code (matched against yang_panel_Sm5_dSv5.png and
# fig. 3's own colourbar labels): temperature = magma over [0,1];
# salinity = Blues over [0, Sm+dSv/2].
YCM_T, YCM_S, YS_LIM = "magma", "Blues", (0.0, 1.5)

FPS = 15


def load_run(d: Path):
    """Load every frame of a run, sorted by time."""
    frames = []
    for p in d.glob("Data_*.h5"):
        m = DATA_RE.search(p.name)
        if not m:
            continue
        with h5py.File(p, "r") as h:
            nx = len(h["grid/xc"]) - 1
            ny = len(h["grid/yc"]) - 1
            nz = len(h["grid/zc"]) - 1
            sl = (slice(0, nz), slice(0, ny), slice(0, nx))
            th = np.asarray(h["Conc/0"][sl], float).mean(axis=0)
            sa = np.asarray(h["Conc/1"][sl], float).mean(axis=0)
            cl = np.asarray(h["VOF/C_L"][sl], float).mean(axis=0)
            x = np.asarray(h["grid/xc"][:nx], float)
            y = np.asarray(h["grid/yc"][:ny], float)
            t = float(np.asarray(h["time"]).flat[0])
        frames.append((t, th, sa, cl, x, y))
    frames.sort(key=lambda f: f[0])
    return frames


def make_video(run: str, out_name: str):
    d = SCR / run
    frames = load_run(d)
    t0, th0, sa0, cl0, x, y = frames[0]
    ext = [x[0], x[-1], y[0], y[-1]]

    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(8.0, 4.0))
    for ax in (ax_t, ax_s):
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.02)

    cmap_t = plt.get_cmap(YCM_T).copy(); cmap_t.set_bad("white")
    cmap_s = plt.get_cmap(YCM_S).copy(); cmap_s.set_bad("white")

    def masked(field, cl):
        return np.ma.array(field, mask=(cl < 0.5))

    im_t = ax_t.imshow(masked(th0, cl0), origin="lower", extent=ext,
                        cmap=cmap_t, vmin=0.0, vmax=1.0, aspect="equal",
                        interpolation="bilinear")
    im_s = ax_s.imshow(masked(sa0, cl0), origin="lower", extent=ext,
                        cmap=cmap_s, vmin=YS_LIM[0], vmax=YS_LIM[1], aspect="equal",
                        interpolation="bilinear")
    cont_t = [ax_t.contour(x, y, cl0, levels=[0.5], colors="k", linewidths=0.8)]
    cont_s = [ax_s.contour(x, y, cl0, levels=[0.5], colors="k", linewidths=0.8)]

    def update(i):
        t, th, sa, cl, _, _ = frames[i]
        im_t.set_data(masked(th, cl))
        im_s.set_data(masked(sa, cl))
        for c in cont_t + cont_s:
            c.remove()
        cont_t[0] = ax_t.contour(x, y, cl, levels=[0.5], colors="k", linewidths=0.8)
        cont_s[0] = ax_s.contour(x, y, cl, levels=[0.5], colors="k", linewidths=0.8)
        return im_t, im_s

    anim = animation.FuncAnimation(fig, update, frames=len(frames), blit=False)
    writer = animation.FFMpegWriter(
        fps=FPS, bitrate=4000, codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    out = VID / out_name
    anim.save(out, writer=writer, dpi=200)
    plt.close(fig)
    print(f"  {out_name}   {len(frames)} frames, t={frames[0][0]:.0f}..{frames[-1][0]:.0f}")


if __name__ == "__main__":
    make_video("YangSalty_07282026/N288", "yang_fields_N288.mp4")
    make_video("YangSalty_07282026/N432", "yang_fields_N432.mp4")
