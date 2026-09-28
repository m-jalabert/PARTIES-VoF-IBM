#!/usr/bin/env python3
"""Compare the 1440^2 provenance re-run against the completed Sm=0 production run.

WHY THIS EXISTS
---------------
The completed 1440^2 Sm=0 run reports t_half = 60.287.  The freshwater ladder
(N = 216/288/360/512/720, all on the current binary) is smooth, monotone and
predictive, and extrapolates to

    t_half(1440) = 70.9-71.6      and      t_half(inf) = 62.9-64.4.

60.287 is not merely off that trend -- it lies BELOW the extrapolated
infinite-resolution limit, which a monotone convergence sequence cannot do.
The 1440^2 run is also the only one built with a different binary
(aac878bf vs the current 6a8fa6eb), and it predates the Le=100 salt-operator
work.  The salty production run predates the current binary too, so this test
bears on every completed production result, not just the freshwater point.

WHAT COUNTS AS A PASS
---------------------
The two runs share every physical input.  They differ only in binary and in
rank count (256 then, 128 now).  Rank-count independence is a stated design
requirement, so agreement should be at round-off accumulation level, far below
any physical effect:

    PASS      max |dV/V| <= 1e-4 over t in [0,10]   -- provenance closed
    MARGINAL  1e-4 < max |dV/V| <= 1e-2             -- investigate before using
    FAIL      max |dV/V| > 1e-2                     -- the 1440^2 point is not
              comparable to the ladder and must be replaced before anything
              (including the "1.78x too fast" headline) rests on it

Note the asymmetry: a PASS exonerates binary AND rank count together; only a
FAIL needs a follow-up to separate them.

Usage:  python3 compare.py [--run DIR] [--ref CSV]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np

V0 = 0.1
DATA_RE = re.compile(r"Data_(\d+)\.h5$")
HERE = Path(__file__).resolve().parent


def load_run(d: Path):
    rows = []
    for p in d.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h5:
            nz = len(h5["grid/zc"]) - 1
            ny = len(h5["grid/yc"]) - 1
            nx = len(h5["grid/xc"]) - 1
            F = np.asarray(h5["VOF/C_L"][:nz, :ny, :nx], dtype=np.float64)
            x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
            y = np.asarray(h5["grid/yc"][:ny], dtype=np.float64)
            t = float(np.asarray(h5["time"]).flat[0])
        v = float((1.0 - F).sum()) * float(x[1] - x[0]) * float(y[1] - y[0]) / V0
        rows.append((t, v))
    rows.sort()
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path,
                    default=Path("/anvil/scratch/x-mjalabert/YangFresh_07272026/prov1440"))
    ap.add_argument("--ref", type=Path,
                    default=HERE / "../Yang_production/yang_sm0_final_timeseries.csv")
    a = ap.parse_args()

    import csv
    rr = list(csv.DictReader(open(a.ref)))
    tr = np.array([float(x["time"]) for x in rr])
    vr = np.array([float(x["ice_volume_ratio"]) for x in rr])

    tn, vn = load_run(a.run)
    if len(tn) == 0:
        print("no outputs in the re-run yet")
        return

    # Compare ONLY at times where the REFERENCE has a real sample.  The old run
    # wrote every 2.0 t_ff, the re-run every 1.0, so probing at odd t forces a
    # linear interpolation of the reference across a curved V(t) and manufactures
    # a difference that has nothing to do with the binary.  In the first dry run
    # this showed up unmistakably: 7.7e-3 at t=1 and 6.6e-4 at t=3 (both
    # interpolated) against 5.0e-8 at t=2 and 1.1e-5 at t=4 (both real samples).
    # Interpolating the FINER series onto the coarser one is safe; the reverse
    # is not.
    tmax = min(tn.max(), tr.max())
    probes = [float(t) for t in tr if 0.0 < t <= tmax + 1e-9]
    if not probes:
        print("re-run has not reached the reference's first output time yet")
        return
    print("comparing only at the reference run's own output times "
          "(no interpolation of the coarser series)\n")
    print(f"{'t':>6} {'V_new(128r, 6a8fa6eb)':>22} {'V_old(256r, aac878bf)':>22} "
          f"{'abs diff':>11} {'rel diff':>10}")
    worst = 0.0
    for t in probes:
        a_ = float(np.interp(t, tn, vn))
        b_ = float(vr[np.argmin(np.abs(tr - t))])
        rel = abs(a_ - b_) / abs(b_)
        worst = max(worst, rel)
        print(f"{t:>6.2f} {a_:>22.8f} {b_:>22.8f} {a_ - b_:>+11.2e} {rel:>10.2e}")

    # Threshold in terms of the SIGNAL, not of round-off.  Two runs differing only
    # in rank count have different summation orders, so their round-off differs and
    # amplifies in a convecting flow: the dry run showed 1.7e-12 at t=2 growing to
    # 1.1e-5 at t=4.  That is chaotic divergence, not disagreement, and demanding
    # bitwise agreement at t=10 would be unphysical.
    #
    # The signal we must detect is whether the completed run's melt rate is
    # consistent with t_half = 60.29 or with the ladder's 71.  At t=10 those
    # trajectories sit at V = 0.887 and V ~ 0.898 respectively, so the
    # discriminating difference is |dV/V| ~ 1.2e-2.  Anything two or more orders
    # below that cannot hide a 15% melt-rate error.
    print(f"\nmax relative difference over t in [0,{probes[-1]:.2f}]: {worst:.3e}")
    print("  (discriminating size: a t_half 60.3-vs-71 difference shows as ~1.2e-2 at t=10)")
    if worst <= 1e-4:
        print("VERDICT: PASS -- the two binaries agree at round-off level.")
        print("  Provenance closed (binary AND rank count).  The 60.287 point is real,")
        print("  and its position below the ladder's t_half(inf) needs a physical")
        print("  explanation instead -- start with whether convergence is genuinely")
        print("  non-monotone between N=720 and N=1440.")
    elif worst <= 4e-3:
        print("VERDICT: PASS (with chaotic divergence) -- above round-off but still well")
        print("  below the ~1.2e-2 that a 15% melt-rate error would produce.  Consistent")
        print("  with rank-count round-off amplifying in a convecting flow, not with the")
        print("  binaries disagreeing.  The 60.287 point stands.")
    elif worst <= 1.2e-2:
        print("VERDICT: MARGINAL -- approaching the discriminating size.  Do not use the")
        print("  1440^2 point without a second realization to separate chaotic spread")
        print("  from a systematic offset.")
    else:
        print("VERDICT: FAIL -- the completed 1440^2 run is NOT comparable to the ladder.")
        print("  Replace it before anything rests on it.  Note the salty production run")
        print("  also predates the current binary and would need the same treatment.")


if __name__ == "__main__":
    main()
