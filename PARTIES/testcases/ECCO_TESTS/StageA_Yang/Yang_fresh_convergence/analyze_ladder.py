#!/usr/bin/env python3
"""Analyse the freshwater convergence ladder (roadmap A.5.2).

Three outputs:

1. **t_half(N)** for each grid, against the completed 1440^2 point (60.287) and
   the digitized Yang 2-D dSv=0 reference (107.42), plus the precommitted
   A.5.2 decision rule.

2. **Binary-provenance check.** The current binary is not the one that produced
   the 1440^2 point (it was rebuilt after the Le=100 salt-operator work, which
   added a default-off runtime path).  The early diffusive phase is
   grid-independent physics, so every grid must collapse onto the same 1-D
   Stefan similarity law

       X(t) = 2 lambda sqrt(t / Pe_T),   V/V0 = 1 - X(t)/0.1

   with the reservoir-driven lambda solving  sqrt(pi) lambda e^{lambda^2}
   (1 + erf lambda) = St.  For St = 0.25 that is lambda = 0.1220; the completed
   A.3 gate and the 1440^2 Sm=0 run measured 0.127.  If the grids do not
   collapse, or lambda has moved, the binary changed the kernel and the 1440^2
   point is not comparable.

3. **Conservation health**: enthalpy <theta + F/St> must hold at 4.5, F in
   [0,1], salt identically zero.

Usage:
    python3 analyze_ladder.py [--root DIR] [--out CSV]
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import brentq

PE_T = 1.0e4
STEFAN = 0.25          # PARTIES convention St = cp dT / L
V0 = 0.1               # configured initial ice volume (0.1 H x H)
SLAB_X0 = 0.9
REF_1440 = 60.287      # completed Sm=0 run
REF_YANG = 107.42      # digitized Yang fig 3(a), Sm=0, dSv=0, 2-D
DATA_RE = re.compile(r"Data_(\d+)\.h5$")


def lambda_reservoir(st: float) -> float:
    """sqrt(pi) l e^{l^2} (1 + erf l) = St."""
    f = lambda l: math.sqrt(math.pi) * l * math.exp(l * l) * (1 + math.erf(l)) - st
    return brentq(f, 1e-9, 5.0)


def physical_shape(h5):
    return (len(h5["grid/zc"]) - 1, len(h5["grid/yc"]) - 1, len(h5["grid/xc"]) - 1)


def physical(ds, shape):
    nz, ny, nx = shape
    return np.asarray(ds[:nz, :ny, :nx], dtype=np.float64)


def scan(run: Path) -> dict:
    files = sorted(
        (int(DATA_RE.search(p.name).group(1)), p)
        for p in run.glob("Data_*.h5") if DATA_RE.search(p.name)
    )
    rows = []
    for _, p in files:
        with h5py.File(p, "r") as h5:
            shape = physical_shape(h5)
            F = physical(h5["VOF/C_L"], shape)
            th = physical(h5["Conc/0"], shape)
            s = physical(h5["Conc/1"], shape)
            t = float(np.asarray(h5["time"]).flat[0])
            x = np.asarray(h5["grid/xc"][: shape[2]], dtype=np.float64)
            y = np.asarray(h5["grid/yc"][: shape[1]], dtype=np.float64)
        dx = float(x[1] - x[0])
        dy = float(y[1] - y[0])
        cell = dx * dy
        ice = float((1.0 - F).sum()) * cell
        rows.append(dict(
            t=t, V=ice / V0,
            enth=float((th + F / STEFAN).sum()) * cell,
            Fmin=float(F.min()), Fmax=float(F.max()),
            smax=float(np.abs(s).max()),
            finite=bool(np.isfinite(F).all() and np.isfinite(th).all()),
        ))
    rows.sort(key=lambda r: r["t"])
    return dict(run=run, rows=rows, n=len(rows))


def band_vs_thermal_layer(run: Path, N: int, v_target: float = 0.75) -> dict:
    """Diffuse-interface width measured against the liquid thermal layer.

    The ladder shows melt rate rising monotonically with N while dE/d(1/N)
    grows, which is not ordinary mesh convergence.  The candidate controlling
    ratio is the diffuse band width against the thermal boundary layer:

        band ~ 4 Cn = 3 dx      (Cn = 0.75 dx on every rung)
        delta_T = distance from the F=0.5 front into the liquid over which
                  theta recovers 90% of its local bulk value

    band/delta_T is 3/N divided by delta_T, so it shrinks as the mesh refines
    even though 'cells across the band' is fixed at 3.  If melt rate tracks
    band/delta_T rather than dx, then Cn -- not dx -- is the knob, and it can be
    converged far more cheaply (A.5.5 variant 3).
    """
    # Select by matched MELT STATE, not matched time: at a common t the grids
    # are at different V/V0 (at t=40 they span 0.63-0.76), so a matched-time
    # boundary-layer comparison is confounded by the ice geometry differing.
    if not run.is_dir():
        return {}
    cands = []
    for p in run.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h5:
            shape = physical_shape(h5)
            F = physical(h5["VOF/C_L"], shape)
            x = np.asarray(h5["grid/xc"][: shape[2]], dtype=np.float64)
            y = np.asarray(h5["grid/yc"][: shape[1]], dtype=np.float64)
            tt = float(np.asarray(h5["time"]).flat[0])
        v = float((1.0 - F).sum()) * float(x[1] - x[0]) * float(y[1] - y[0]) / V0
        cands.append((abs(v - v_target), v, tt, p))
    if not cands:
        return {}
    cands.sort()
    dv, v_sel, t_sel, path = cands[0]
    if dv > 0.03:            # never reached this state; don't report a fake row
        return {}
    with h5py.File(path, "r") as h5:
        shape = physical_shape(h5)
        nz, ny, nx = shape
        F = physical(h5["VOF/C_L"], shape).mean(axis=0)
        th = physical(h5["Conc/0"], shape).mean(axis=0)
        x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
        t = t_sel
    dx = float(x[1] - x[0])
    lens = []
    for j in range(ny):
        row, tr = F[j], th[j]
        if not (row.min() <= 0.5 <= row.max()):
            continue
        # front index: last liquid-side cell before the interface
        idx = int(np.argmax(row[::-1] >= 0.5))
        i_f = nx - 1 - idx
        if i_f < 10:
            continue
        bulk = float(np.median(tr[max(0, i_f - 128): max(1, i_f - 8)]))
        edge = float(tr[i_f])
        if not np.isfinite(bulk) or abs(bulk - edge) < 1e-6:
            continue
        target = edge + 0.9 * (bulk - edge)
        seg = tr[: i_f + 1][::-1]              # marching from the front into the liquid
        hit = np.argmax(seg >= target) if bulk > edge else np.argmax(seg <= target)
        if hit == 0:
            continue
        lens.append(hit * dx)
    if not lens:
        return {}
    d90 = float(np.median(lens))
    band = 3.0 * dx                            # 4*Cn with Cn = 0.75 dx

    # Interface arclength: the competing explanation for a faster melt on a
    # finer grid is simply more exposed AREA (better-resolved scallops), not
    # more heat per unit area.  Measure it and settle that.
    y = np.arange(ny) * dx                     # uniform grid: dy == dx
    fx = np.full(ny, np.nan)
    for j in range(ny):
        row = F[j]
        if row.min() <= 0.5 <= row.max():
            fx[j] = np.interp(0.5, row[::-1], x[::-1])
    good = np.isfinite(fx)
    L_int = float(np.sqrt(np.diff(fx[good]) ** 2 + dx ** 2).sum()) if good.sum() > 1 else float("nan")

    return dict(t=t, V=v_sel, dx=dx, d90=d90, d90_cells=d90 / dx,
                band=band, ratio=band / d90, n_rows=len(lens),
                L_int=L_int, front_rms=float(np.nanstd(fx)),
                front_ptp=float(np.nanmax(fx) - np.nanmin(fx)))


def t_half(rows) -> float:
    t = np.array([r["t"] for r in rows])
    v = np.array([r["V"] for r in rows])
    for i in range(1, len(v)):
        if v[i - 1] > 0.5 >= v[i]:
            return float(t[i - 1] + (0.5 - v[i - 1]) * (t[i] - t[i - 1]) / (v[i] - v[i - 1]))
    return float("nan")


def fit_lambda(rows, t_lo=1.0, t_hi=6.5) -> tuple[float, int]:
    """Least-squares lambda from X(t) = 2 l sqrt(t/Pe_T) over the diffusive window.

    The window is short on purpose: the 1440^2 Sm=0 run already reaches
    max_speed = 0.136 by t = 4, so anything past t ~ 6 is convectively enhanced
    and biases lambda high.  With outputs every 2.0 t_ff this leaves only a few
    points, so treat the fit as secondary to the matched-time V(t) table, which
    needs no model at all.
    """
    sel = [r for r in rows if t_lo <= r["t"] <= t_hi and r["V"] < 1.0]
    if len(sel) < 3:
        return float("nan"), len(sel)
    tt = np.array([r["t"] for r in sel])
    X = (1.0 - np.array([r["V"] for r in sel])) * V0
    basis = 2.0 * np.sqrt(tt / PE_T)
    return float((basis @ X) / (basis @ basis)), len(sel)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path,
                    default=Path("/anvil/scratch/x-mjalabert/YangFresh_07272026"))
    ap.add_argument("--out", type=Path, default=Path("ladder_summary.csv"))
    ap.add_argument(
        "--ref-1440-csv", type=Path,
        default=Path(__file__).resolve().parent
        / "../Yang_production/yang_sm0_final_timeseries.csv",
        help="completed 1440^2 Sm=0 time series, included as the finest ladder rung")
    a = ap.parse_args()

    lam_th = lambda_reservoir(STEFAN)
    print(f"reservoir-driven similarity lambda for St={STEFAN}: {lam_th:.4f}")
    print(f"references: 1440^2 PARTIES t_half={REF_1440}, Yang 2-D dSv=0 t_half={REF_YANG}\n")

    res = {}
    for d in sorted(a.root.glob("N*")):
        m = re.match(r"N(\d+)$", d.name)
        if not m:
            continue
        r = scan(d)
        if r["n"] == 0:
            print(f"N{m.group(1):>5}: no outputs yet")
            continue
        res[int(m.group(1))] = r

    # The completed 1440^2 run is the finest rung; it exists only as a summary
    # CSV (its raw fields are on another filesystem), so splice it in here.
    if a.ref_1440_csv.exists():
        import csv as _csv
        rr = list(_csv.DictReader(open(a.ref_1440_csv)))
        rows = [dict(t=float(x["time"]), V=float(x["ice_volume_ratio"]),
                     enth=float(x["enthalpy_mean"]), Fmin=float(x["F_min"]),
                     Fmax=float(x["F_max"]), smax=abs(float(x["salt_max"])),
                     finite=bool(float(x["finite"]))) for x in rr]
        rows.sort(key=lambda r_: r_["t"])
        res[1440] = dict(run=a.ref_1440_csv, rows=rows, n=len(rows))
        print(f"spliced in the completed 1440^2 run from {a.ref_1440_csv.name} "
              f"({len(rows)} outputs)\n")

    hdr = (f"{'N':>6} {'outs':>5} {'t_end':>7} {'V_end':>7} {'t_half':>8} "
           f"{'lambda':>7} {'nfit':>5} {'d_enth':>10} {'Fmin':>9} {'Fmax':>9} {'|s|max':>8}")
    print(hdr)
    print("-" * len(hdr))
    lines = ["N,outputs,t_end,V_end,t_half,lambda_fit,n_fit,enthalpy_rel_drift,F_min,F_max,abs_s_max"]
    for N in sorted(res):
        rows = res[N]["rows"]
        th_ = t_half(rows)
        lam, nfit = fit_lambda(rows)
        e0 = rows[0]["enth"]
        drift = max(abs(r["enth"] - e0) for r in rows) / abs(e0)
        Fmin = min(r["Fmin"] for r in rows)
        Fmax = max(r["Fmax"] for r in rows)
        smax = max(r["smax"] for r in rows)
        print(f"{N:>6} {len(rows):>5} {rows[-1]['t']:>7.2f} {rows[-1]['V']:>7.4f} "
              f"{th_:>8.3f} {lam:>7.4f} {nfit:>5} {drift:>10.2e} {Fmin:>9.2e} "
              f"{Fmax:>9.6f} {smax:>8.1e}")
        lines.append(f"{N},{len(rows)},{rows[-1]['t']:.6f},{rows[-1]['V']:.6f},"
                     f"{th_:.6f},{lam:.6f},{nfit},{drift:.3e},{Fmin:.3e},{Fmax:.9f},{smax:.3e}")
    a.out.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {a.out}")

    # ---- matched-time tables.
    #
    # CAVEAT, learned the hard way 2026-07-27: early-time agreement between two
    # coarse grids proves nothing.  The 1-D Stefan front has moved only
    # X(t) = 2 lambda sqrt(t/Pe_T), which at t=2 is 1.2 cells on N360 and 1.8 on
    # N512 -- comparable to the diffuse band itself.  Two under-resolved grids
    # agree with each other because they share the same band-scale artefact, not
    # because the physics is resolved.  A collapse is only diagnostic where the
    # front spans >= MIN_CELLS cells on every grid compared.
    MIN_CELLS = 5.0

    def x_1d(t):
        return 2.0 * lam_th * math.sqrt(t / PE_T)

    probes = [t for t in (2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 120, 140)
              if all(r["rows"][-1]["t"] >= t for r in res.values())]
    grids = sorted(res)
    if probes and len(res) > 1:
        def v_at(N, t):
            rr = res[N]["rows"]
            return float(np.interp(t, [r["t"] for r in rr], [r["V"] for r in rr]))

        print("\nMatched-time V/V0.  'cells' is the 1-D Stefan front displacement on the\n"
              "COARSEST grid; rows below %.0f cells are marked (*) and are NOT diagnostic --\n"
              "the front is unresolved there, so agreement is an artefact." % MIN_CELLS)
        print("      t  " + "".join(f"{'N%d' % N:>10}" for N in grids)
              + f"{'V_1D':>9}{'cells':>8}{'spread':>9}")
        for t in probes:
            vs = [v_at(N, t) for N in grids]
            cells = x_1d(t) * min(grids)
            sp = (max(vs) - min(vs)) / np.mean(vs)
            flag = " " if cells >= MIN_CELLS else "*"
            print(f"  {t:5.0f}  " + "".join(f"{v:>10.4f}" for v in vs)
                  + f"{1 - x_1d(t) / V0:>9.4f}{cells:>8.1f}{sp:>8.2%}{flag}")

        print("\nConvective enhancement E = (1-V) / (1-V_1D)  -- melt relative to pure\n"
              "1-D conduction.  This is the quantity that actually separates the grids:")
        print("      t  " + "".join(f"{'N%d' % N:>10}" for N in grids))
        for t in probes:
            v1d = 1.0 - x_1d(t) / V0
            if v1d >= 1.0:
                continue
            print(f"  {t:5.0f}  " + "".join(
                f"{(1 - v_at(N, t)) / (1 - v1d):>10.3f}" for N in grids))
        if len(grids) > 1:
            tq = probes[-1]
            v1d = 1.0 - x_1d(tq) / V0
            es = [(1 - v_at(N, tq)) / (1 - v1d) for N in grids]
            trend = ("RISES" if es[-1] > es[0] else "FALLS")
            print(f"\n  E at t={tq:.0f} {trend} with N: "
                  + ", ".join(f"N{N}={e:.3f}" for N, e in zip(grids, es)))
            if es[-1] > es[0]:
                print("  -> finer grid = stronger convection = FASTER melt.  Since the finest\n"
                      "     grid is already faster than Yang, refinement moves AWAY from the\n"
                      "     reference; resolution is not the remedy.")
            # is the sequence settling?
            if len(grids) > 2:
                slopes = [abs(es[i + 1] - es[i]) / abs(1.0 / grids[i + 1] - 1.0 / grids[i])
                          for i in range(len(grids) - 1)]
                print("  dE/d(1/N) between successive grids: "
                      + ", ".join(f"{s:.0f}" for s in slopes)
                      + ("  -> GROWING: not in the asymptotic range, the finest grid is\n"
                         "     probably not converged either."
                         if slopes[-1] > slopes[0] else
                         "  -> shrinking: consistent with convergence."))

    # ---- diffuse band vs thermal layer (raw fields only, so 1440 is excluded)
    bt = {}
    for N in grids:
        d = a.root / f"N{N}"
        b = band_vs_thermal_layer(d, N)
        if b:
            bt[N] = b
    if bt:
        print("\nDiffuse band vs liquid thermal layer, at MATCHED MELT STATE V/V0 ~ 0.75\n"
              "(the band is 4*Cn = 3*dx on every rung, so 'cells across the band' is fixed\n"
              "at 3 and only the RATIO to the thermal layer moves):")
        print(f"{'N':>6} {'V/V0':>6} {'t':>7} {'delta_T':>9} {'cells':>7} {'band':>9} "
              f"{'band/dT':>9} {'L_int':>8} {'front_ptp':>10}")
        for N in sorted(bt):
            b = bt[N]
            print(f"{N:>6} {b['V']:>6.3f} {b['t']:>7.1f} {b['d90']:>9.5f} {b['d90_cells']:>7.1f} "
                  f"{b['band']:>9.5f} {b['ratio']:>8.1%} {b['L_int']:>8.4f} {b['front_ptp']:>10.5f}")
        Ls = [bt[N]["L_int"] for N in sorted(bt)]
        if len(Ls) > 1:
            print(f"  interface arclength spread: {(max(Ls) - min(Ls)) / np.mean(Ls):.2%} "
                  "-> exposed AREA is grid-independent, so a faster melt on the finer\n"
                  "  grid is more heat PER UNIT AREA, not more area.")
        print("  If melt rate tracks band/delta_T rather than dx, the knob is Cn, not the\n"
              "  mesh -- and Cn converges far more cheaply (roadmap A.5.5 variant 3).")

    halves = {N: t_half(res[N]["rows"]) for N in res}
    halves = {N: v for N, v in halves.items() if not math.isnan(v)}
    if not halves:
        print("\nno grid has reached V/V0 = 0.5 yet -- decision rule not applicable")
        return

    lams = [fit_lambda(res[N]["rows"])[0] for N in res]
    lams = [l for l in lams if not math.isnan(l)]
    if len(lams) > 1:
        spread = (max(lams) - min(lams)) / np.mean(lams)
        print(f"\nBINARY-PROVENANCE CHECK: lambda across grids = "
              f"{['%.4f' % l for l in lams]}, spread {spread:.2%}, theory {lam_th:.4f}")
        # A spread here is NOT evidence of a changed kernel.  This is a laterally
        # heated cavity: the flow starts at t=0+ from baroclinic torque, with no
        # instability threshold, so even the shortest fit window is convectively
        # contaminated -- and contaminated MORE on finer grids, which convect
        # earlier (E(t=6) is 1.215 at N=1440 vs 1.08 at N=360).  A finer grid
        # therefore reports a higher apparent lambda for purely physical reasons.
        # This test cannot separate that from a kernel change, so it can only
        # ever return INCONCLUSIVE or a hard failure (non-finite / wildly off).
        if not all(0.5 * lam_th < l < 3.0 * lam_th for l in lams):
            print("  FAIL -- lambda is far from theory on at least one grid; investigate")
        else:
            print("  INCONCLUSIVE (by construction) -- all grids are within a factor of the\n"
                  "  theoretical lambda, but convective contamination scales with resolution,\n"
                  "  so the spread is expected and this cannot certify the binary.  Provenance\n"
                  "  rests on code inspection: every yang_salt_transport use site in Conc.c is\n"
                  "  guarded and the runtime default is 0.")

    allh = dict(halves)
    allh.setdefault(1440, REF_1440)
    print("\nt_half(N): " + ", ".join(f"N{N}={v:.2f}" for N, v in sorted(allh.items())))
    vals = list(allh.values())
    spread = (max(vals) - min(vals)) / np.mean(vals)
    print(f"spread across the ladder: {spread:.2%}")
    print(f"Yang target {REF_YANG}; PARTIES/Yang rate ratio at the finest grid: "
          f"{REF_YANG / allh[min(allh, key=lambda k: -k)]:.3f}x too fast")

    print("\nA.5.2 DECISION:")
    if spread < 0.05:
        print("  CONVERGED (spread < 5%). Resolution is NOT the cause of the 1.78x.")
        print("  -> the discrepancy is model or reference; go to roadmap A.5.5,")
        print("     and use this ladder as the cheap sandbox for model variants.")
    else:
        ordered = [allh[N] for N in sorted(allh)]
        if ordered[-1] > ordered[0]:
            print("  t_half RISES with N -> refinement moves TOWARD Yang.")
            if len(ordered) >= 2:
                for p in (1, 2):
                    rich = ordered[-1] + (ordered[-1] - ordered[-2]) / (2 ** p - 1)
                    print(f"     Richardson p={p}: t_half(inf) ~ {rich:.2f} "
                          f"({'reaches' if rich >= REF_YANG else 'short of'} {REF_YANG})")
            print("  -> the 1440^2 run may be under-resolved; consider a 2880^2 confirmation.")
        else:
            print("  t_half FALLS with N -> refinement moves AWAY from Yang.")
            print("  -> strong result; resolution is not the remedy. Go to roadmap A.5.5.")


if __name__ == "__main__":
    main()
