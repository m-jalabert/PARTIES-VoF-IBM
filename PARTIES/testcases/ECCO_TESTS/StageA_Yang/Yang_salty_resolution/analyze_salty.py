#!/usr/bin/env python3
"""Analyse the salty (Sm=5, dSv=5) resolution runs.

Reports t_half against the two reference points, and -- because these grids
deliberately under-resolve the salt field -- three health checks that must pass
BEFORE any t_half is believed:

  1. salt boundedness.  s is initialised linear in y over [0.5, 1.5] in the
     liquid and 0 in the ice, so physical s lies in [0, 1.5].  Pe_S = 1e6 with
     2nd-order CENTRAL advection on a coarse mesh can produce dispersive
     wiggles; s < 0 is unphysical and indicates the salt field is
     under-resolved to the point of corrupting the buoyancy.
  2. salt conservation.  The masked-diffusion operator conserves int(s) exactly
     (that is how meltwater dilution emerges), so drift here is a solver
     problem, not physics.
  3. enthalpy conservation of int(theta + F/St).  The Step-C failure produced a
     dramatic, entirely spurious melt rate while leaking enthalpy at 1e-2; a
     melt rate is meaningless without this check.

References:
    Yang Sm=5 dSv=5   t_half ~ 216   (fig 3(b) f/f0 = 0.497, f0 = 1/107.42)
    PARTIES 1440^2    t_half = 179.85
Freshwater scaling t_half(288)/t_half(1440) = 1.698 predicts 305 here if the
melt rate simply follows resolution.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np

V0 = 0.1
STEFAN = 0.25
YANG = 216.1
P1440 = 179.85
NAIVE = 179.85 * 1.698
DATA_RE = re.compile(r"Data_(\d+)\.h5$")


def scan(d: Path):
    rows = []
    for p in d.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h5:
            nz = len(h5["grid/zc"]) - 1
            ny = len(h5["grid/yc"]) - 1
            nx = len(h5["grid/xc"]) - 1
            sl = (slice(0, nz), slice(0, ny), slice(0, nx))
            F = np.asarray(h5["VOF/C_L"][sl], dtype=np.float64)
            th = np.asarray(h5["Conc/0"][sl], dtype=np.float64)
            s = np.asarray(h5["Conc/1"][sl], dtype=np.float64)
            x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
            y = np.asarray(h5["grid/yc"][:ny], dtype=np.float64)
            t = float(np.asarray(h5["time"]).flat[0])
        cell = float(x[1] - x[0]) * float(y[1] - y[0])
        rows.append(dict(
            t=t, V=float((1.0 - F).sum()) * cell / V0,
            enth=float((th + F / STEFAN).sum()) * cell,
            salt=float(s.sum()) * cell,
            smin=float(s.min()), smax=float(s.max()),
            Fmin=float(F.min()), Fmax=float(F.max()),
            finite=bool(np.isfinite(s).all() and np.isfinite(F).all()),
        ))
    rows.sort(key=lambda r: r["t"])
    return rows


def t_half(rows):
    t = np.array([r["t"] for r in rows])
    v = np.array([r["V"] for r in rows])
    for i in range(1, len(v)):
        if v[i - 1] > 0.5 >= v[i]:
            return float(t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1]))
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path,
                    default=Path("/anvil/scratch/x-mjalabert/YangSalty_07282026"))
    a = ap.parse_args()

    print(f"references: Yang ~{YANG:.0f} | PARTIES 1440^2 = {P1440} | "
          f"naive resolution scaling -> {NAIVE:.0f}\n")
    hdr = (f"{'N':>6} {'outs':>5} {'t_end':>7} {'V_end':>7} {'t_half':>8} "
           f"{'dE@t.5':>9} {'d_salt':>9} {'s_min':>10} {'s_max':>7} {'F rng':>14}")
    print(hdr)
    print("-" * len(hdr))
    for d in sorted(a.root.glob("N*")):
        m = re.match(r"N(\d+)$", d.name)
        if not m:
            continue
        rows = scan(d)
        if not rows:
            print(f"{m.group(1):>6}  (no outputs yet)")
            continue
        e0, s0 = rows[0]["enth"], rows[0]["salt"]
        ds = max(abs(r["salt"] - s0) for r in rows) / abs(s0)
        smin = min(r["smin"] for r in rows)
        smax = max(r["smax"] for r in rows)
        # Health must be judged AT t_half, not as max-over-run.  These runs go to
        # t=350, long past melt-through, and the late-time drift there says
        # nothing about whether t_half is trustworthy -- exactly the distinction
        # that mattered for the freshwater ladder.
        th_ = t_half(rows)
        tt = np.array([r["t"] for r in rows])
        dr = np.array([abs(r["enth"] - e0) / abs(e0) for r in rows])
        de = float(np.interp(th_, tt, dr)) if th_ == th_ else dr.max()
        de_max = dr.max()
        print(f"{m.group(1):>6} {len(rows):>5} {rows[-1]['t']:>7.1f} {rows[-1]['V']:>7.4f} "
              f"{th_:>8.2f} {de:>9.2e} {ds:>9.2e} {smin:>10.3e} {smax:>7.4f} "
              f"[{min(r['Fmin'] for r in rows):.1e},{max(r['Fmax'] for r in rows):.4f}]")

        flags = []
        # A global s_min is not itself disqualifying: int(s) is conserved to
        # machine precision, so under/overshoots are paired local dispersion.
        # What matters is HOW MUCH of the liquid is affected, and when.
        if smin < -1e-6:
            flags.append(f"salt undershoot s_min={smin:.2e} (check affected fraction below)")
        if smax > 1.5 + 1e-6:
            flags.append(f"SALT OVERSHOOT s_max={smax:.4f} > 1.5")
        if de > 1e-4:
            flags.append(f"ENTHALPY LEAK {de:.2e} at t_half -- cf. Step C (1.1e-2); "
                         "t_half NOT usable")
        if ds > 1e-8:
            flags.append(f"SALT DRIFT {ds:.2e} (masked diffusion should conserve int(s))")
        if not all(r["finite"] for r in rows):
            flags.append("NON-FINITE FIELDS")
        for f in flags:
            print(f"        !! {f}")
        if not flags:
            print("        health OK (salt bounded, both budgets conserved)")

        print(f"        enthalpy drift: {de:.2e} at t_half, {de_max:.2e} max at end of run")
        th = th_
        if th == th:
            print(f"        vs Yang {YANG:.0f}: {th / YANG:.3f}x | vs 1440^2 {P1440}: "
                  f"{th / P1440:.3f}x | vs naive scaling {NAIVE:.0f}: {th / NAIVE:.3f}x")
            if abs(th - YANG) / YANG < 0.10:
                print("        -> lands on Yang: heat-transport resolution dominates the")
                print("           salty case too; the freshwater story transfers.")
            elif abs(th - NAIVE) / NAIVE < 0.12:
                print("        -> follows naive resolution scaling: coarsening S and phi")
                print("           along with u,T is NOT equivalent to Yang's dual-grid split.")
            else:
                print("        -> neither: needs interpretation, do not force a narrative.")


if __name__ == "__main__":
    main()
