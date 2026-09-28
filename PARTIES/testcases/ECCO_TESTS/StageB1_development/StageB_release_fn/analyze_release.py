#!/usr/bin/env python3
"""
B.1.3 functional release test -- analyser.

The default-off regression proved the release code is a no-op when disabled.
This asks the opposite question, and it is the one that matters for ECCO:
with the mechanism ON, does the grain

  (1) stay EXACTLY still while frozen in the ice,
  (2) come free when its surroundings have actually melted (not before, not
      never), and
  (3) hand off to the flow WITHOUT AN IMPULSE?

(3) is the subtle one.  The grain accumulates a hydrodynamic force history
while held, so an unramped release would dump that accumulated force in as a
single kick.  The ramp is supposed to prevent that; this checks whether it
does, by comparing the acceleration across release against the free-settling
acceleration scale the grain should be feeling.

Usage:  python3 analyze_release.py [rundir]
"""
import sys, os, glob, re
import numpy as np

RUN = sys.argv[1] if len(sys.argv) > 1 else "run64"

# --- physical scales, READ FROM THE DECK -------------------------------------
# Parsed rather than hardcoded so this analyses any release deck, not just the
# one it was written against.


def _deck(rundir):
    """Scalar inputs from parties.inp (last assignment wins; comments stripped)."""
    v = {}
    path = os.path.join(rundir, "parties.inp")
    if not os.path.exists(path):
        path = os.path.join(os.path.dirname(rundir.rstrip("/")) or ".",
                            "parties.inp")
    with open(path) as f:
        for ln in f:
            ln = ln.split("#")[0].split("//")[0].strip()
            if "=" not in ln or ln.startswith("["):
                continue
            k, _, val = ln.partition("=")
            v[k.strip()] = val.strip()
    return v


def _f(v, k, default=None):
    try:
        return float(v[k])
    except Exception:
        return default


_V = _deck(RUN)
with open(os.path.join(
        RUN if os.path.exists(os.path.join(RUN, "p_mobile.inp"))
        else os.path.dirname(RUN.rstrip("/")) or ".", "p_mobile.inp")) as _f_:
    _lines = [l for l in _f_.read().split("\n") if l.strip()]
_x0, _y0, _z0, _R = (float(t) for t in _lines[1].split()[:4])

R        = _R
D        = 2 * R
RHO_S    = _f(_V, "rho_s", 2.0)
G        = 1.0
Y0       = _y0
YMIN     = _f(_V, "ymin", 0.0)
YMAX     = _f(_V, "ymax", 1.0)
Y_INT0   = YMIN + _f(_V, "vof_slab_x0", 0.5) * (YMAX - YMIN)
F_REL    = _f(_V, "F_release", 0.9)
RAMP     = int(_f(_V, "release_ramp_steps", 10))
# Free-settling acceleration of the grain the instant it is released, in the
# absence of any hydrodynamic response: a = g*(rho_s - rho_f)/rho_s.
A_FREE   = G * (RHO_S - 1.0) / RHO_S

print(f"deck: R={R}, grain y0={Y0}, rho_s={RHO_S}, y_interface0={Y_INT0}, "
      f"F_release={F_REL}, ramp={RAMP}")


def load_release(path):
    t, ID, phi, trel, ramp, y, u0, u1, u2 = [], [], [], [], [], [], [], [], []
    with open(path) as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.strip().split(",")
            if len(p) < 9:
                continue
            t.append(float(p[0])); ID.append(int(p[1])); phi.append(float(p[2]))
            trel.append(float(p[3])); ramp.append(int(p[4])); y.append(float(p[5]))
            u0.append(float(p[6])); u1.append(float(p[7])); u2.append(float(p[8]))
    return {k: np.asarray(v) for k, v in
            dict(t=t, ID=ID, phi=phi, trel=trel, ramp=ramp,
                 y=y, u0=u0, u1=u1, u2=u2).items()}


def melt_front(rundir):
    """Ice/ocean interface height vs time, from the F=0.5 level of the
    horizontally averaged phase field -- away from the grain's column so the
    grain itself does not bias it."""
    import h5py
    out = []
    files = sorted(glob.glob(os.path.join(rundir, "Data_*.h5")),
                   key=lambda p: int(re.search(r"_(\d+)\.h5$", p).group(1)))
    for fn in files:
        with h5py.File(fn) as f:
            t = float(f["time"][0])
            F = f["VOF/F"][0]          # (ny, nx)
            yc = f["grid/yc"][:]
            nx = F.shape[1]
            # exclude the middle third: that is the grain's column
            cols = np.r_[0:nx // 3, 2 * nx // 3:nx]
            prof = F[:, cols].mean(axis=1)
            # F = 1 (water) below, 0 (ice) above -> first crossing of 0.5
            k = np.argmax(prof < 0.5)
            if k == 0:
                out.append((t, np.nan)); continue
            f0, f1 = prof[k - 1], prof[k]
            yi = yc[k - 1] + (0.5 - f0) / (f1 - f0) * (yc[k] - yc[k - 1])
            out.append((t, yi))
    return np.array(out)


def main():
    rel_path = os.path.join(RUN, "release.dat")
    if not os.path.exists(rel_path):
        print(f"FAIL: {rel_path} not found -- release diagnostic never wrote.")
        print("      Either F_release <= 0 (mechanism off) or the build lacks")
        print("      VOF_IBM + LAG_PARTICLE_RESOLVED.")
        return 1

    d = load_release(rel_path)
    n = len(d["t"])
    print(f"release.dat: {n} steps, t in [{d['t'][0]:.4g}, {d['t'][-1]:.4g}]\n")

    checks = []

    # ---- (1) did it fire, and when? -----------------------------------------
    fired = d["trel"] > 0
    if not fired.any():
        print("RESULT: the grain was NEVER released.")
        print(f"  max phi_liq reached = {d['phi'].max():.4f} "
              f"(threshold {F_REL})")
        print("  -> either the run is too short, or the criterion cannot reach")
        print("     the threshold, which would be a real defect.")
        checks.append(("release fires", False))
        t_rel = None
    else:
        t_rel = d["trel"][fired][0]
        i_rel = int(np.argmax(fired))
        print(f"[1] RELEASED at t = {t_rel:.6f}  (step index {i_rel})")
        print(f"    phi_liq at release = {d['phi'][i_rel]:.4f}  "
              f"(threshold {F_REL})")
        checks.append(("release fires", True))

    # ---- (2) exactly locked before release ----------------------------------
    if t_rel is not None:
        pre = d["t"] < t_rel
        umax = max(np.abs(d["u0"][pre]).max(), np.abs(d["u1"][pre]).max(),
                   np.abs(d["u2"][pre]).max()) if pre.any() else 0.0
        dy = abs(d["y"][pre][-1] - Y0) if pre.any() else 0.0
        ok = (umax == 0.0) and (dy < 1e-12)
        print(f"\n[2] LOCKED BEFORE RELEASE:  max|U| = {umax:.3e}, "
              f"|dy| = {dy:.3e}   -> {'PASS' if ok else 'FAIL'}")
        print("    (must be EXACTLY zero: the lock zeroes U and U_old, so any")
        print("     nonzero value means the lock is leaking)")
        checks.append(("locked before release", ok))

    # ---- (3) impulse-free handoff -------------------------------------------
    if t_rel is not None:
        post = d["t"] >= t_rel
        tp, up = d["t"][post], d["u1"][post]
        if len(tp) > 3:
            dt = np.diff(tp)
            acc = np.diff(up) / np.where(dt > 0, dt, np.nan)
            acc = acc[np.isfinite(acc)]
            a_peak = np.abs(acc).max() if len(acc) else 0.0
            ratio = a_peak / A_FREE
            # A grain released into quiescent fluid cannot accelerate faster
            # than free settling.  Allow 2x for the added-mass/pressure
            # transient; anything beyond that is an impulse.
            ok = ratio < 2.0
            print(f"\n[3] IMPULSE-FREE HANDOFF:")
            print(f"    peak |dU_y/dt| after release = {a_peak:.4g}")
            print(f"    free-settling scale a_free   = {A_FREE:.4g}")
            print(f"    ratio = {ratio:.2f}x   -> {'PASS' if ok else 'FAIL'} "
                  f"(bar: < 2x)")
            # velocity at the end of the ramp
            in_ramp = post & (d["ramp"] > 0)
            if in_ramp.any():
                print(f"    ramp active for {in_ramp.sum()} steps "
                      f"(release_ramp_steps = {RAMP})")
                print(f"    U_y at ramp end = {d['u1'][in_ramp][-1]:.4e}")
            checks.append(("impulse-free handoff", ok))

    # ---- (4) does it actually settle? ---------------------------------------
    if t_rel is not None:
        post = d["t"] >= t_rel
        if post.sum() > 10:
            dy_fall = d["y"][post][-1] - d["y"][post][0]
            uy_mean = d["u1"][post][len(d["u1"][post]) // 2:].mean()
            ok = dy_fall < -1e-3 and uy_mean < 0
            print(f"\n[4] SETTLES AFTER RELEASE:")
            print(f"    dy over the post-release window = {dy_fall:+.4f}")
            print(f"    mean U_y (second half)          = {uy_mean:+.4e}")
            print(f"    -> {'PASS' if ok else 'FAIL'} (grain must fall)")
            checks.append(("settles after release", ok))

    # ---- (5) is the criterion physically sensible? --------------------------
    print(f"\n[5] CRITERION TRAJECTORY:")
    print(f"    phi_liq: {d['phi'][0]:.4f} -> {d['phi'].max():.4f}")
    # monotone-ish: fraction of steps where it decreases materially
    dphi = np.diff(d["phi"])
    frac_down = (dphi < -1e-6).sum() / max(len(dphi), 1)
    print(f"    fraction of steps where phi_liq decreases = {frac_down:.3f}")
    print("    (melting is one-way here, so phi_liq should be near-monotone;")
    print("     a large decreasing fraction would mean the shell integral is")
    print("     picking up advected liquid rather than the melt state)")

    # cross-check against the melt front
    try:
        mf = melt_front(RUN)
        good = mf[np.isfinite(mf[:, 1])]
        if len(good):
            print(f"\n[6] MELT FRONT CROSS-CHECK "
                  f"(F=0.5, away from the grain's column):")
            print(f"    y_interface: {good[0,1]:.4f} (t={good[0,0]:.2f}) "
                  f"-> {good[-1,1]:.4f} (t={good[-1,0]:.2f})")
            print(f"    initial interface (deck)  = {Y_INT0:.4f}")
            print(f"    grain bottom              = {Y0 - R:.4f}")
            if t_rel is not None:
                yi_rel = np.interp(t_rel, good[:, 0], good[:, 1])
                # Ice sits ABOVE the interface and the front RISES as the ice
                # melts, so exposure grows with y_int.  Measure it from the
                # bottom of the grain in grain diameters:
                #   0 = front just reaching the grain's underside
                #   1 = front level with the grain's top (fully uncovered)
                #  >1 = front has passed the grain entirely
                expo = (yi_rel - (Y0 - R)) / D
                print(f"    interface at release      = {yi_rel:.4f}")
                print(f"    exposure at release       = {expo:.2f} grain "
                      f"diameters above the grain's underside")
                if expo >= 1.0:
                    print(f"    -> the front had ALREADY cleared the whole grain by "
                          f"{expo - 1.0:.2f} d when it fired:")
                    print("       F_release = %.2f is a LATE (conservative) trigger."
                          % F_REL)
                elif expo <= 0.0:
                    print("    -> fired while still fully buried: too early.")
                else:
                    print(f"    -> fired with the grain {100*(1-expo):.0f} % still "
                          "covered by ice.")
    except Exception as e:
        print(f"\n[6] melt-front cross-check unavailable: {e}")

    # ---- (7) meltwater tracer conservation (B.1.4) --------------------------
    # Conc_add_meltwater_RHS injects the tracer at exactly melt_src, the same
    # field that pays the latent heat.  With no sink and no-flux walls that
    # forces  d/dt int(C_mw) = d/dt int(F), i.e. the tracer total must track
    # the ice volume lost cell for cell.  Any drift is a transport error --
    # which is precisely what has to be ruled out before entrainment numbers
    # from this tracer mean anything.
    try:
        import h5py
        files = sorted(glob.glob(os.path.join(RUN, "Data_*.h5")),
                       key=lambda p_: int(re.search(r"_(\d+)\.h5$", p_).group(1)))
        rows = []
        for fn in files:
            with h5py.File(fn) as f:
                if "Conc/2" not in f:
                    rows = []
                    break
                tt = float(f["time"][0])
                rows.append((tt, f["Conc/2"][0].sum(), f["VOF/F"][0].sum()))
        if len(rows) > 2:
            a = np.array(rows)
            mw = a[:, 1] - a[0, 1]          # tracer gained
            ice = a[:, 2] - a[0, 2]         # liquid gained = ice lost
            print("\n[7] MELTWATER TRACER CONSERVATION (B.1.4):")
            print(f"    {'t':>7} {'int C_mw':>13} {'int dF':>13} {'rel drift':>11}")
            for tt, m, ic in zip(a[:, 0], mw, ice):
                rel = abs(m - ic) / abs(ic) if abs(ic) > 1e-12 else 0.0
                print(f"    {tt:7.2f} {m:13.5e} {ic:13.5e} {rel:11.3e}")
            fin = abs(mw[-1] - ice[-1]) / max(abs(ice[-1]), 1e-30)
            ok = fin < 1e-3
            print(f"    final relative drift = {fin:.3e}  -> "
                  f"{'PASS' if ok else 'FAIL'} (bar: < 1e-3)")
            checks.append(("meltwater conservation", ok))
        else:
            print("\n[7] meltwater tracer not present in this run "
                  "(NConc < 3 or meltwater_tracer = 0) -- skipped")
    except Exception as e:
        print(f"\n[7] meltwater conservation check unavailable: {e}")

    # ---- verdict ------------------------------------------------------------
    print("\n" + "=" * 62)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    allok = all(ok for _, ok in checks) and len(checks) > 0
    print(f"\nVERDICT: {'PASS' if allok else 'FAIL'}"
          " -- interface-triggered release is functional"
          if allok else "\nVERDICT: FAIL")
    print("=" * 62)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
