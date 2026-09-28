#!/usr/bin/env python3
"""Re-fit the Roquet quadratic EOS to polar (ECCO) conditions.  Roadmap B.1.1.

The Stage-A eos_* coefficients were tuned to Yang's warm, nearly fresh lab
window (0-20 C, Sm = 5).  ECCO is a different corner of the (T,S) plane --
near-freezing water at S ~ 30-35 g/kg -- and they do not transfer.

Fits exactly the form the solver evaluates in Velocity.c:

    b(theta, s) = -betaT * |theta - Tmd0 - Tmd_slope*s|^q + betaS * s

TWO THINGS MAKE THE NAIVE FIT WRONG, and both are handled here.

(1) The temperature scale must be the superheat above the LOCAL freezing
    point, not the gap to the freshwater melting point.  At S = 34 the
    freezing point is -1.85 C, so measuring theta from 0 C gives a
    temperature scale of ~0.15 K instead of the physical ~2 K and the whole
    nondimensionalization collapses.  Here theta = (T - Tf(S_ocean))/dT with
    dT the prescribed superheat, so theta = 0 is the ice-side temperature and
    theta = 1 the far-field ocean -- the convention the ice initial conditions
    already use.  The solver's own T_melt input then carries the offset to the
    s = 0 melting point, and is reported below.

(2) The parabola's vertex is NOT determined by data in the polar box.  At
    S >= 28 the temperature of maximum density lies BELOW the freezing point
    (see the table at the bottom), so within the liquid range the quadratic is
    evaluated on one monotonic flank, far from its vertex -- where it is
    locally indistinguishable from a straight line.  Fitting betaT, Tmd0 and
    Tmd_slope freely is therefore ill-posed: they trade off along a valley and
    a least-squares solver happily returns betaS ~ 1e8 with a 1e-10 prefactor
    and a 0.2% residual.  That fit is numerically excellent and physically
    meaningless.

    The fix is to PIN the vertex to the true TEOS-10 density-maximum locus,
    which is a property of seawater and is known exactly, and fit only the two
    buoyancy amplitudes betaT and betaS -- which the data DO determine.
"""
from __future__ import annotations

import argparse

import gsw
import numpy as np


def t_freeze(SA, p=0.0):
    """Conservative-temperature freezing point (degC)."""
    return float(gsw.CT_freezing(SA, p, 0.0))


def t_maxdensity(SA, p=0.0):
    """Conservative temperature of maximum density (degC), by refined search."""
    T = np.linspace(-15.0, 10.0, 2501)
    rho = gsw.rho(np.full_like(T, SA), T, p)
    i = int(np.argmax(rho))
    lo, hi = T[max(i - 1, 0)], T[min(i + 1, len(T) - 1)]
    T2 = np.linspace(lo, hi, 2001)
    rho2 = gsw.rho(np.full_like(T2, SA), T2, p)
    return float(T2[int(np.argmax(rho2))])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--S-ocean", type=float, default=34.0)
    ap.add_argument("--S-min", type=float, default=28.0)
    ap.add_argument("--S-max", type=float, default=35.0)
    ap.add_argument("--superheat", type=float, default=2.0,
                    help="far-field temperature above the local freezing point, K")
    ap.add_argument("--p", type=float, default=0.0, help="pressure, dbar")
    a = ap.parse_args()

    p, S_ocean = a.p, a.S_ocean
    Tf_ocean = t_freeze(S_ocean, p)
    T_ref = Tf_ocean                      # theta = 0
    dT = a.superheat                      # theta = 1 at the far field
    dS = S_ocean                          # s = 1 at the far field
    T_ocean = T_ref + dT

    print("=" * 74)
    print("ECCO polar EOS re-fit  (roadmap B.1.1)")
    print("=" * 74)
    print(f"box:    S = {a.S_min}..{a.S_max} g/kg, "
          f"T = Tf(S) .. Tf(S) + {dT} K,  p = {p} dbar")
    print(f"scales: theta = (T - {T_ref:.4f} C)/{dT} K,   s = S/{dS} g/kg")
    print(f"        theta = 0 -> ice-side ({T_ref:.4f} C), "
          f"theta = 1 -> ocean ({T_ocean:.4f} C)")

    # --- sample the box ----------------------------------------------------
    SA_ax = np.linspace(a.S_min, a.S_max, 80)
    SAg, CTg = [], []
    for s_ in SA_ax:
        tf = t_freeze(s_, p)
        for t_ in np.linspace(tf, tf + dT, 80):
            SAg.append(s_)
            CTg.append(t_)
    SAg, CTg = np.array(SAg), np.array(CTg)
    rho = gsw.rho(SAg, CTg, p)

    theta = (CTg - T_ref) / dT
    s = SAg / dS

    # --- (2) pin the vertex to the true density-maximum locus --------------
    Smd = np.linspace(a.S_min, a.S_max, 25)
    Tmd = np.array([t_maxdensity(x, p) for x in Smd])
    b_md, a_md = np.polyfit(Smd, Tmd, 1)      # T_md(S) = a_md + b_md*S
    md_res = float(np.max(np.abs(np.polyval([b_md, a_md], Smd) - Tmd)))

    Tmd0 = (a_md - T_ref) / dT
    slope = b_md * dS / dT

    print(f"\ndensity-maximum locus from TEOS-10: "
          f"T_md(S) = {a_md:.4f} + {b_md:.6f}*S  (max dev {md_res:.4f} K)")

    # --- model comparison --------------------------------------------------
    #
    # A FREE ADDITIVE CONSTANT MUST BE ALLOWED.  Buoyancy enters the momentum
    # equation as a force, so a spatially uniform offset is absorbed into the
    # pressure gradient and is dynamically irrelevant.  Omitting it makes every
    # model look catastrophically bad (rms/std ~ 99%) purely because the mean
    # is unconstrained -- that is a fitting artefact, not physics.
    rho_ref = float(gsw.rho(S_ocean, T_ocean, p))
    b_target = -(rho - rho_ref) / rho_ref     # buoyancy anomaly, dimensionless
    sd = float(np.std(b_target))
    one = np.ones_like(s)
    dth = theta - Tmd0 - slope * s

    def fit_model(A, name):
        c, *_ = np.linalg.lstsq(A, b_target, rcond=None)
        r = A @ c - b_target
        rms = float(np.sqrt(np.mean(r ** 2)))
        print(f"  {name:<42} rms/std = {rms/sd:8.4%}   "
              f"cond = {np.linalg.cond(A):8.3g}")
        return c, rms / sd

    print("\n--- model comparison (free constant allowed) ---")
    c_lin, e_lin = fit_model(np.column_stack([one, theta, s]),
                             "LINEAR   c0 + c1*theta + c2*s")
    c_roq, e_roq = fit_model(np.column_stack([one, -(dth * dth), s]),
                             "ROQUET   quadratic, vertex pinned")
    _, e_full = fit_model(
        np.column_stack([one, theta, s, theta ** 2, theta * s, s ** 2]),
        "full quadratic (attainable ceiling)")

    betaT_d, betaS_d = float(c_roq[1]), float(c_roq[2])

    print("\n--- RECOMMENDATION ---")
    print(f"A LINEAR EOS is sufficient for this box: {e_lin:.3%} of the")
    print(f"buoyancy variance is unexplained, against {e_roq:.3%} for the")
    print("quadratic.  The quadratic buys a factor ~2.5 on an already")
    print("negligible residual, and it buys it dishonestly: the fitted")
    print(f"betaT is {betaT_d:+.4e}" +
          (" -- NEGATIVE, which the Roquet form does not admit.\n"
           "The parabola is being used far from its vertex as a proxy for a\n"
           "linear trend, which is exactly the misuse the vertex lies outside\n"
           "the liquid range to warn about."
           if betaT_d < 0 else "."))
    print("\nSo for ECCO: do NOT compile EOS_NONLINEAR.  Use the linear")
    print("buoyancy path with these coefficients:")
    print(f"    d(b)/d(theta) = {float(c_lin[1]):+.6e}")
    print(f"    d(b)/d(s)     = {float(c_lin[2]):+.6e}")
    print(f"    haline/thermal buoyancy ratio over the box = "
          f"{abs(float(c_lin[2])/float(c_lin[1])):.0f}")
    print("  (signs are in the b = -(rho-rho_ref)/rho_ref convention; check")
    print("   them against the momentum assembly before entering richardson.)")
    print("\nFor reference, had the quadratic been used, normalized to")
    print("betaT = 1 as in Stage A:")
    if betaT_d != 0:
        print(f"    eos_betaS = {betaS_d/betaT_d:.6f}, "
              f"eos_Tmd0 = {Tmd0:.6f}, eos_Tmd_slope = {slope:.6f}")

    # --- the solver's other phase-change inputs ----------------------------
    T_melt_nd = (0.0 - T_ref) / dT
    tf_lo, tf_hi = t_freeze(a.S_min, p), t_freeze(a.S_max, p)
    dTf_dS = (tf_hi - tf_lo) / (a.S_max - a.S_min)
    liq_nd = -dTf_dS * dS / dT

    print("\n--- matching [phase_change] inputs ---")
    print(f"T_melt         = {T_melt_nd:.6f}   "
          f"(freshwater melting point 0 C in these units)")
    print(f"liquidus_slope = {liq_nd:.6f}   (= -dTf/dS * dS/dT, "
          f"dTf/dS = {dTf_dS:.6f} K/(g/kg))")
    print(f"  Stage A used 0.014.  The polar value is ~{liq_nd/0.014:.0f}x larger:")
    print(f"  the freezing-point depression across the full salinity range")
    print(f"  ({-dTf_dS*dS:.2f} K) now dwarfs the {dT} K superheat, whereas in")
    print(f"  Yang's warm case the superheat was ~20 K.  Salinity therefore")
    print(f"  controls melting far more strongly here than in the benchmark.")

    # --- the qualitative point --------------------------------------------
    print("\n--- density maximum vs freezing point ---")
    print(f"{'S (g/kg)':>10} {'T_freeze':>10} {'T_maxdens':>11}  location")
    for s_ in (0.0, 5.0, 10.0, 20.0, 24.0, 28.0, 34.0, 35.0):
        tf = t_freeze(s_, p) if s_ > 0 else 0.0
        tmd = t_maxdensity(max(s_, 1e-6), p)
        print(f"{s_:>10.1f} {tf:>10.4f} {tmd:>11.4f}  "
              f"{'INSIDE liquid range' if tmd > tf else 'below freezing'}")
    print("\nAt polar salinity the density maximum lies BELOW the freezing point,")
    print("so liquid seawater has no density maximum: density rises")
    print("monotonically as T falls.  Yang's Sm = 5 case sits on the other side")
    print("of that transition -- which is precisely why the benchmark shows the")
    print("non-monotonic melt-rate behaviour it does.")
    print("\nConsequence for ECCO: the anomalous-convection regime validated in")
    print("Stage A does NOT carry over; meltwater is unambiguously buoyant here.")
    print("What Stage A validated is the machinery -- Stefan coupling, salt")
    print("transport, penalization -- not the regime.  State it that way.")


if __name__ == "__main__":
    main()
