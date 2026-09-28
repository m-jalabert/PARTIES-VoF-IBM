#!/usr/bin/env python3
"""Digitize Yang et al. (2023) figure 3 from the arXiv LaTeX source figure.

Source: arXiv:2302.02357 e-print bundle, ``fig3.pdf``
        (sha256 a6bdd059f3fc986ad9b9916e46c163aff6141d6e081731b3b3d9dd7de2c30bd9),
        archived at docs/papers/yang2023_arxiv_source/fig3.pdf.
        JFM 969 R2 figure 4 is the same figure.

Panel (a): V(t)/V0 vs t for Sm = 0, 5, 10, 15 g/kg.
Panel (b): fbar/fbar_0 vs Sm, colour-coded by dSv; circles 2-D, squares 3-D.

Because the figure is vector art rendered by us at 900 dpi, the only error is
line width, not screen capture.  The V = 0.5 dashed line lands on the frame
mid-height to 1e-4, which is the built-in calibration check.

Usage:  python3 digitize_yang_fig3.py [path/to/fig3.pdf] [outdir]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

DPI = 900
HERE = Path(__file__).resolve().parent
PDF = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    HERE / "../../../../../docs/papers/yang2023_arxiv_source/fig3.pdf").resolve()
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE


def render(pdf: Path) -> np.ndarray:
    stem = OUT / "_fig3_render"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-png", "-singlefile",
                    str(pdf), str(stem)], check=True)
    return np.asarray(Image.open(str(stem) + ".png").convert("RGB")).astype(int)


def group(idx, tol=4):
    out, cur = [], [idx[0]]
    for v in idx[1:]:
        if v - cur[-1] <= tol:
            cur.append(v)
        else:
            out.append(int(np.mean(cur)))
            cur = [v]
    out.append(int(np.mean(cur)))
    return out


def frame(sub):
    dark = sub.sum(axis=2) < 200
    h, w = dark.shape
    vs = group(np.where(dark.sum(axis=0) > 0.5 * h)[0])
    hs = group(np.where(dark.sum(axis=1) > 0.4 * w)[0])
    return vs, hs, dark


# ---------------------------------------------------------------- panel (a)
def panel_a(a):
    W = a.shape[1]
    sub = a[:, : int(W * 0.50)]
    vs, hs, _ = frame(sub)
    x0, x1, y0, y1 = vs[0], vs[-1], hs[0], hs[-1]
    px2t = lambda p: (p - x0) / (x1 - x0) * 200.0
    py2v = lambda p: 1.0 + (p - y0) / (y1 - y0) * (0.0 - 1.0)
    assert abs(py2v((y0 + y1) / 2) - 0.5) < 1e-3, "panel (a) calibration failed"

    nominal = {"Sm=0": (31, 59, 99), "Sm=5": (168, 207, 232),
               "Sm=10": (244, 177, 131), "Sm=15": (192, 57, 43)}
    pix = sub.reshape(-1, 3)
    colours = {}
    for name, c0 in nominal.items():
        sel = pix[np.abs(pix - np.array(c0)).sum(axis=1) < 90]
        colours[name] = np.median(sel, axis=0).astype(int)

    lx = (x0, x0 + int(0.62 * (x1 - x0)))          # legend box interior
    ly = (y0 + int(0.55 * (y1 - y0)), y1)

    out = {}
    for name, c in colours.items():
        ts, vv = [], []
        for px in range(x0 + 3, x1 - 2):
            col = sub[y0 + 2: y1 - 2, px]
            m = np.where(np.abs(col - c).sum(axis=1) < 60)[0]
            if not len(m):
                continue
            pys = m + y0 + 2
            if lx[0] <= px <= lx[1]:
                pys = pys[(pys < ly[0]) | (pys > ly[1])]
                if not len(pys):
                    continue
            ts.append(px2t(px))
            vv.append(py2v(np.median(pys)))
        out[name] = (np.array(ts), np.array(vv))
    return out


def metrics(t, v):
    o = np.argsort(t)
    t, v = t[o], v[o]
    th = np.nan
    for i in range(1, len(v)):
        if v[i - 1] > 0.5 >= v[i]:
            th = t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1])
            break
    return th, t, v


# ---------------------------------------------------------------- panel (b)
def panel_b(a):
    W = a.shape[1]
    sub = a[:, int(W * 0.50):]
    vs, hs, dark = frame(sub)
    x0, x1, y0, y1 = vs[0], vs[1], hs[0], hs[-1]

    rows = sorted(r for r in group(np.where(dark[:, x0 + 3: x0 + 20].sum(axis=1) > 10)[0])
                  if y0 + 15 < r < y1 + 5)
    cols = sorted(c for c in group(np.where(dark[y1 - 20: y1 - 3, :].sum(axis=0) > 10)[0])
                  if x0 - 5 < c < x1 - 15)
    ycoef = np.polyfit(rows, [0.2 + 0.2 * i for i in range(len(rows))][::-1], 1)
    xcoef = np.polyfit(cols, [5.0 * i for i in range(len(cols))], 1)

    interior = np.zeros(sub.shape[:2], bool)
    interior[y0 + 3: y1 - 3, x0 + 3: x1 - 3] = True
    cand = interior & ~(sub.sum(axis=2) > 720) & ~(sub.sum(axis=2) < 260)
    lab, _ = ndimage.label(cand)

    pts = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        ys, xs = sl
        hh, ww = ys.stop - ys.start, xs.stop - xs.start
        m = lab[sl] == i
        if not (18 <= hh <= 90 and 18 <= ww <= 90):
            continue
        if m.sum() < 0.55 * hh * ww or abs(hh - ww) > 0.30 * max(hh, ww):
            continue
        c = np.median(sub[sl][m], axis=0)
        pts.append((xcoef[0] * (xs.start + xs.stop) / 2 + xcoef[1],
                    ycoef[0] * (ys.start + ys.stop) / 2 + ycoef[1], c))

    red = lambda c: c[0] > 150 and c[1] < 90 and c[2] < 90
    minima = sorted((s, f) for s, f, c in pts if red(c))
    blues = [(s, f, c) for s, f, c in pts if not red(c)]

    arr = np.array([c for _, _, c in blues])
    used, series = np.zeros(len(blues), bool), []
    for i in np.argsort(arr.sum(axis=1))[::-1]:            # lightest (small dSv) first
        if used[i]:
            continue
        sel = (np.abs(arr - arr[i]).sum(axis=1) < 45) & ~used
        used |= sel
        series.append(sorted((blues[j][0], blues[j][1]) for j in np.where(sel)[0]))
    return series, minima


def main():
    a = render(PDF)
    curves = panel_a(a)

    print("=== panel (a): V(t)/V0, 2-D, dSv = 0 ===")
    rows, half = [], {}
    for name in ("Sm=0", "Sm=5", "Sm=10", "Sm=15"):
        th, t, v = metrics(*curves[name])
        half[name] = th
        print(f"  {name:6s} t_half = {th:7.2f}   V(200) = {np.interp(199.0, t, v):.4f}")
        for tq in range(0, 201, 5):
            rows.append((name.split("=")[1], tq, float(np.interp(tq, t, v))))
    print("  f_bar ratios vs Sm=0: " + ", ".join(
        f"{n}={half['Sm=0'] / half[n]:.4f}" for n in ("Sm=5", "Sm=10", "Sm=15")))

    p = OUT / "yang_fig3a_reference.csv"
    with open(p, "w") as fh:
        fh.write("# Yang et al. JFM 969 R2 (2023) / arXiv:2302.02357 figure 3(a) [= JFM fig 4(a)]\n")
        fh.write("# 2-D simulations, dSv = 0 (identified by matching the f_bar ratios to the\n")
        fh.write("# lightest circle series of panel (b); see YANG_REFERENCE_CORRECTION_REPORT.md)\n")
        fh.write("# digitized from the arXiv vector figure at 900 dpi; uncertainty ~0.005 in V/V0\n")
        fh.write("Sm_g_per_kg,time_t_ff,ice_volume_ratio\n")
        for sm, tq, v in rows:
            fh.write(f"{sm},{tq},{v:.5f}\n")
    print(f"  wrote {p}")

    series, minima = panel_b(a)
    print("\n=== panel (b): f_bar/f_bar_0, lightest (dSv=0) first ===")
    p = OUT / "yang_fig3b_reference.csv"
    with open(p, "w") as fh:
        fh.write("# Yang figure 3(b) [= JFM fig 4(b)]: normalized melt rate f_bar/f_bar_0.\n")
        fh.write("# series_rank 0 = lightest = dSv 0; ranks increase with dSv.\n")
        fh.write("# The digitizer recovers the well-separated series; overlapping markers\n")
        fh.write("# in the two darkest series are partially missed.  2-D circles only.\n")
        fh.write("series_rank,Sm_g_per_kg,f_ratio\n")
        for k, s in enumerate(series):
            print(f"  rank {k}: " + "  ".join(f"({sm:g},{f:.3f})" for sm, f in s))
            for sm, f in s:
                fh.write(f"{k},{sm:.2f},{f:.4f}\n")
    print(f"  wrote {p}")
    print("\n  red minimum markers:", [(round(s, 2), round(f, 3)) for s, f in minima])


if __name__ == "__main__":
    main()
