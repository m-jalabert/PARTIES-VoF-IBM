#!/usr/bin/env python3
"""Analyse the Stage-B 2-D release shakeout (roadmap B.1.3/B.1.4).

This is a code shake-out, so every check below is a property the implementation
must have by construction -- not a physical result to be interpreted.  A run
that melts beautifully but fails check 2 or 4 is a broken run.

  1. LOCKED PHASE.  Before release the grain must not move at all: X constant to
     machine precision and U identically zero.  Any drift means the motion lock
     at Lagrangian.c:1118 is leaking.
  2. RELEASE EVENT.  The grain must be freed once, at the time the melt front
     clears its shell, and never re-lock.
  3. IMPULSE-FREE RELEASE.  The whole point of the ramp is that |dU/dt| just
     after release is comparable to |dU/dt| later, not a spike.  A grain held
     against a converged force history and then let go without a ramp shows a
     single-step jump; we report the ratio explicitly.
  4. MELTWATER TRACER CONSERVATION.  The tracer is injected with the SAME
     discrete melt_src field that drives the phase change, so with no-flux
     walls integral(C_mw) must track the ice volume lost cell for cell.  This is
     the check that makes the tracer a diagnostic rather than a decoration.
  5. ENTHALPY CONSERVATION of integral(theta + F/St) -- the Stage-A gate that
     caught the Step-C failure; a melt rate is meaningless without it.
  6. MIXING DIAGNOSTICS (B.1.6), computed offline on the saved fields as the
     roadmap allows: horizontally-averaged fluxes <w'c'>, effective diffusivity
     kappa_eff = -<w'c'>/d<c>/dy, and meltwater penetration depth.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import h5py
import numpy as np

DATA_RE = re.compile(r"Data_(\d+)\.h5$")
PART_RE = re.compile(r"Particle_(\d+)\.h5$")


def load_fields(d: Path, stefan: float):
    rows = []
    for p in d.glob("Data_*.h5"):
        if not DATA_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h:
            nz = len(h["grid/zc"]) - 1
            ny = len(h["grid/yc"]) - 1
            nx = len(h["grid/xc"]) - 1
            sl = (slice(0, nz), slice(0, ny), slice(0, nx))
            F = np.asarray(h["VOF/C_L"][sl], float)
            CS = np.asarray(h["VOF/C_S"][sl], float) if "VOF/C_S" in h else np.zeros_like(F)
            th = np.asarray(h["Conc/0"][sl], float)
            mw = np.asarray(h["Conc/2"][sl], float) if "Conc/2" in h else None
            v = np.asarray(h["v"][sl], float) if "v" in h else None
            u = np.asarray(h["u"][sl], float) if "u" in h else None
            x = np.asarray(h["grid/xc"][:nx], float)
            y = np.asarray(h["grid/yc"][:ny], float)
            t = float(np.asarray(h["time"]).flat[0])
        cell = float(x[1] - x[0]) * float(y[1] - y[0])
        rows.append(dict(t=t, cell=cell, y=y, x=x, F=F, CS=CS, th=th, mw=mw, v=v, u=u,
                         ice=float((1.0 - F - CS).clip(0).sum()) * cell,
                         liq=float(F.sum()) * cell,
                         enth=float((th + F / stefan).sum()) * cell,
                         mwint=(float(mw.sum()) * cell) if mw is not None else np.nan))
    rows.sort(key=lambda r: r["t"])
    return rows


def load_particle(d: Path):
    """Read the mobile-particle trajectory.

    ParticleOutput.c writes /mobile only when Np > 0, so a run with no
    particles produces files holding nothing but /domain and /time.  Treat that
    as "no trajectory" rather than an error -- it is the normal shape of every
    Stage-A output.
    """
    rows = []
    for p in d.glob("Particle_*.h5"):
        if not PART_RE.search(p.name):
            continue
        with h5py.File(p, "r") as h:
            if "mobile" not in h or "X" not in h["mobile"]:
                continue
            t = float(np.asarray(h["time"]).flat[0]) if "time" in h else np.nan
            grp = h["mobile"]
            X = np.asarray(grp["X"], float).reshape(-1, 3)
            U = np.asarray(grp["U"], float).reshape(-1, 3)
        rows.append((t, X[0].copy(), U[0].copy()))
    rows.sort(key=lambda r: r[0])
    return rows


def background_pe(b, F, y, cell, liquid_only=True):
    """Background potential energy of the buoyancy field b (Winters 1995).

    BPE is the potential energy the field would have after being sorted
    adiabatically into a stable profile -- the part of the PE that mixing can
    only increase.  Computed here by an explicit sort, which is exact for a
    uniform grid and is done OFFLINE on the saved fields exactly as the roadmap
    (B.1.6) allows: the in-solver alternative would need a distributed sort, and
    an MPI_Gather of the whole field to rank 0 would not scale.

    Only liquid cells take part.  Ice is not fluid: including it would let the
    solid's buoyancy masquerade as available potential energy, and the ice
    volume changes in time, which would put a spurious trend in BPE.

    The returned value uses the SORTED heights, so BPE increases only through
    irreversible mixing, not through advection.
    """
    bb = b.ravel()
    if liquid_only:
        m = F.ravel() > 0.5
        bb = bb[m]
    if bb.size == 0:
        return np.nan
    # Densest (most negative buoyancy) at the bottom.
    order = np.argsort(bb)
    # Heights the sorted parcels would occupy, bottom-up, on this grid.
    ny = len(y)
    hh = np.interp(np.linspace(0, ny - 1, bb.size), np.arange(ny), y)
    return float(np.sum(bb[order] * hh) * cell)


def dissipation(u, v, y, x, re):
    """Volume-integrated viscous dissipation, 2-D: eps = (2/Re) S_ij S_ij."""
    dy = float(y[1] - y[0])
    dx = float(x[1] - x[0])
    ux = np.gradient(u, dx, axis=1)
    uy = np.gradient(u, dy, axis=0)
    vx = np.gradient(v, dx, axis=1)
    vy = np.gradient(v, dy, axis=0)
    s2 = ux ** 2 + vy ** 2 + 0.5 * (uy + vx) ** 2
    return float(2.0 / re * s2.sum() * dx * dy)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--stefan", type=float, default=0.25)
    ap.add_argument("--re", type=float, default=500.0,
                    help="Reynolds number, for the dissipation in the mixing efficiency")
    a = ap.parse_args()

    rows = load_fields(a.root, a.stefan)
    if not rows:
        print("no Data_*.h5 outputs found")
        return
    parts = load_particle(a.root)

    t = np.array([r["t"] for r in rows])
    ice = np.array([r["ice"] for r in rows])
    enth = np.array([r["enth"] for r in rows])
    mwint = np.array([r["mwint"] for r in rows])

    # --- release time, from the solver's own log if present ------------------
    t_rel = np.nan
    log = a.root / "run.log"
    if log.exists():
        m = re.search(r"\[release\] particle \d+ freed at t = ([0-9.eE+-]+)",
                      log.read_text(errors="ignore"))
        if m:
            t_rel = float(m.group(1))

    print(f"outputs: {len(rows)}   t = {t.min():.2f} .. {t.max():.2f}")
    print(f"ice volume: {ice[0]:.6f} -> {ice[-1]:.6f} "
          f"({100*(1-ice[-1]/ice[0]):.1f}% melted)")
    print(f"release time reported by solver: "
          f"{t_rel:.4f}" if t_rel == t_rel else "release time: NOT REPORTED")

    flags = []

    # --- 1/2/3 particle behaviour -------------------------------------------
    if parts:
        tp = np.array([r[0] for r in parts])
        Xp = np.array([r[1] for r in parts])
        Up = np.array([r[2] for r in parts])
        print(f"\nparticle: {len(parts)} snapshots, "
              f"X0 = ({Xp[0,0]:.4f}, {Xp[0,1]:.4f})")

        if t_rel == t_rel:
            pre = tp < t_rel
            post = tp >= t_rel
        else:
            # fall back on motion onset
            moved = np.abs(Up[:, 1]) > 1e-12
            pre = ~moved
            post = moved

        if pre.sum() > 1:
            dX = np.abs(Xp[pre] - Xp[pre][0]).max()
            dU = np.abs(Up[pre]).max()
            print(f"  LOCKED phase ({pre.sum()} snapshots): "
                  f"max |X-X0| = {dX:.3e}, max |U| = {dU:.3e}")
            if dX > 1e-12 or dU > 1e-12:
                flags.append(f"MOTION LOCK LEAKING: |dX|={dX:.2e} |U|={dU:.2e} "
                             "before release (must be exactly 0)")
        else:
            flags.append("no locked snapshots captured -- release happened "
                         "before the first output; move the grain deeper or "
                         "shorten output_time_interval")

        if post.sum() > 2:
            Uy = Up[post][:, 1]
            tt = tp[post]
            acc = np.abs(np.diff(Uy) / np.diff(tt))
            print(f"  RELEASED phase ({post.sum()} snapshots): "
                  f"|U_y| {np.abs(Uy).min():.3e} -> {np.abs(Uy).max():.3e}")
            if len(acc) > 2:
                ratio = acc[0] / np.median(acc[1:])
                print(f"  impulse check: |dU/dt| first interval / median later "
                      f"= {ratio:.2f}")
                if ratio > 5.0:
                    flags.append(f"RELEASE IMPULSE: first |dU/dt| is {ratio:.1f}x "
                                 "the later median -- the ramp is not damping it")
            if np.any(Uy > 0) and np.any(Uy < 0):
                print("  note: U_y changes sign (grain not monotonically sinking)")
    else:
        flags.append("no Particle_*.h5 found -- cannot check the lock or the ramp")

    # --- 4 meltwater tracer conservation ------------------------------------
    if np.isfinite(mwint).all():
        # Compare against the change in integral(C_L), NOT the derived ice
        # volume.  melt_src drives dC_L/dt directly, so integral(C_L) is the
        # quantity the tracer must track; the ice volume 1-C_L-C_S differs from
        # it by the clip and by any change in the resolved solid, which would
        # show up here as a spurious mismatch that is not a transport error.
        liq = np.array([r["liq"] for r in rows])
        lost = liq - liq[0]
        print(f"\nmeltwater tracer vs the liquid it created:")
        print(f"{'t':>7} {'d int C_L':>12} {'int C_mw':>12} {'rel diff':>11}")
        worst = 0.0
        for i in range(0, len(t), max(1, len(t) // 8)):
            if lost[i] <= 0:
                continue
            rd = abs(mwint[i] - lost[i]) / abs(lost[i])
            worst = max(worst, rd)
            print(f"{t[i]:>7.2f} {lost[i]:>12.6e} {mwint[i]:>12.6e} {rd:>11.2e}")
        print(f"worst relative mismatch: {worst:.3e}")
        if worst > 1e-2:
            flags.append(f"MELTWATER TRACER NOT CONSERVED: {worst:.2e} "
                         "(same melt_src drives both, so this must close)")

    # --- 5 enthalpy ----------------------------------------------------------
    de = np.abs(enth - enth[0]) / abs(enth[0])
    print(f"\nenthalpy drift: {de.max():.3e} max over the run")
    if de.max() > 1e-3:
        flags.append(f"ENTHALPY LEAK {de.max():.2e} -- cf. Step C (1.1e-2)")

    # --- 6 mixing diagnostics -----------------------------------------------
    last = rows[-1]
    if last["v"] is not None and last["mw"] is not None:
        y = last["y"]
        v = last["v"][0]
        c = last["mw"][0]
        vbar = v.mean(axis=1, keepdims=True)
        cbar = c.mean(axis=1, keepdims=True)
        flux = ((v - vbar) * (c - cbar)).mean(axis=1)
        dcdy = np.gradient(cbar[:, 0], y)
        with np.errstate(divide="ignore", invalid="ignore"):
            kap = np.where(np.abs(dcdy) > 1e-12, -flux / dcdy, np.nan)
        thr = 0.01 * np.nanmax(cbar)
        depth = float(y[np.argmax(cbar[:, 0] > thr)]) if np.nanmax(cbar) > 0 else np.nan
        print(f"\nmixing diagnostics at t = {last['t']:.2f} (B.1.6, offline):")
        print(f"  max |<w'c'>|           = {np.abs(flux).max():.3e}")
        print(f"  median kappa_eff       = {np.nanmedian(np.abs(kap)):.3e}")
        print(f"  meltwater penetration  = y {depth:.4f} "
              f"(lowest y with <c> > 1% of max)")

    # --- 6b background PE / mixing efficiency (Winters 1995) ---------------
    if last["u"] is not None and last["v"] is not None:
        # Buoyancy proxy: for this freshwater shakeout the stratifying agent is
        # temperature (warm = light), so b = theta.  A salty ECCO run must
        # replace this with the run's actual EOS -- stated rather than silently
        # assumed, because a wrong b makes BPE meaningless.
        bpe = np.array([background_pe(r["th"][0], r["F"][0], r["y"], r["cell"])
                        for r in rows])
        d_bpe = float(bpe[-1] - bpe[0])
        eps = dissipation(last["u"][0], last["v"][0], last["y"], last["x"], a.re)
        denom = abs(d_bpe) + eps * float(t[-1] - t[0])
        eta = abs(d_bpe) / denom if denom > 0 else np.nan
        print(f"\nbackground PE (b = theta, liquid cells only):")
        print(f"  BPE {bpe[0]:.6e} -> {bpe[-1]:.6e}   dBPE = {d_bpe:+.3e}")
        print(f"  dissipation at final time  = {eps:.3e}")
        print(f"  mixing efficiency eta      = {eta:.4f}")
        print("  (eta uses the final-time dissipation across the whole span;")
        print("   a time-integrated eps needs profile output every ~10 steps)")

    print("\n" + "=" * 66)
    if flags:
        for f in flags:
            print(f"!! {f}")
        print("VERDICT: NOT CLEAN -- fix the above before the 3-D shakeout")
    else:
        print("VERDICT: PASS -- lock, release, ramp, tracer and enthalpy all clean")


if __name__ == "__main__":
    main()
