#!/usr/bin/env python3
"""
B.1.1 -- re-fit the Roquet quadratic EOS to POLAR conditions against TEOS-10.

The Stage-A eos_* coefficients are Yang's laboratory values (Cb = 0.011
kg/m^3/K^2, b0 = 0.77 kg/m^3/(g/kg), T0 = 4 C, cS = -0.25 K/(g/kg)), tuned to a
warm freshwater tank.  ECCO is the opposite corner of the (T,S) plane: near
freezing, S ~ 28-35 g/kg.  This fits the SAME functional form the code
implements

    rho(T,S) = rho0 - Cb*|T - T0 - cS*S|^q + b0*S        (q = 2)

to TEOS-10 over the polar box, and converts the result to the nondimensional
eos_* inputs.

Code form (Velocity.c:2252):  b = -betaT*|theta - Tmd0 - slope*s|^q + betaS*s

No code change -- only input values.
"""
import numpy as np
import gsw
from scipy.optimize import least_squares

# --- the polar box -----------------------------------------------------------
# Absolute Salinity and Conservative Temperature, TEOS-10's native variables.
SA_LO, SA_HI = 28.0, 35.0        # g/kg
CT_LO, CT_HI = -2.2, 2.0         # deg C
P_REF        = 0.0               # dbar (near-surface); sensitivity at 50 dbar below

N = 60
SA = np.linspace(SA_LO, SA_HI, N)
CT = np.linspace(CT_LO, CT_HI, N)
SAg, CTg = np.meshgrid(SA, CT, indexing="ij")

# Mask points below the freezing line -- they are not ocean.
CTf = gsw.CT_freezing(SAg, P_REF, 0.0)
mask = CTg >= CTf
print(f"polar box: SA in [{SA_LO}, {SA_HI}] g/kg, CT in [{CT_LO}, {CT_HI}] C, "
      f"p = {P_REF} dbar")
print(f"  {mask.sum()} of {mask.size} grid points are above the freezing line\n")

s = SAg[mask]
t = CTg[mask]


def rho_teos(SA_, CT_, p):
    return gsw.rho(SA_, CT_, p)


rho_ref = rho_teos(s, t, P_REF)


# --- model ------------------------------------------------------------------
def rho_model(par, S_, T_, q=2.0):
    rho0, Cb, T0, cS, b0 = par
    return rho0 - Cb * np.abs(T_ - T0 - cS * S_) ** q + b0 * S_


def resid(par):
    return rho_model(par, s, t) - rho_ref


# Yang's lab values as the starting point (and as the "before" row).
YANG = np.array([1000.0, 0.011, 4.0, -0.25, 0.77])

fit = least_squares(resid, YANG, method="lm", max_nfev=200000)
rho0, Cb, T0, cS, b0 = fit.x

r_yang = rho_model(YANG, s, t) - rho_ref
r_fit = fit.fun

print("=" * 74)
print("QUADRATIC ROQUET FIT (q = 2), dimensional")
print("=" * 74)
print(f"{'':22s} {'Yang (lab)':>14s} {'polar re-fit':>14s}")
print(f"{'rho0  [kg/m^3]':22s} {YANG[0]:14.5f} {rho0:14.5f}")
print(f"{'Cb    [kg/m^3/K^2]':22s} {YANG[1]:14.6f} {Cb:14.6f}")
print(f"{'T0    [C]':22s} {YANG[2]:14.4f} {T0:14.4f}")
print(f"{'cS    [K/(g/kg)]':22s} {YANG[3]:14.4f} {cS:14.4f}")
print(f"{'b0    [kg/m^3/(g/kg)]':22s} {YANG[4]:14.4f} {b0:14.4f}")
print()
print(f"{'RMS error [kg/m^3]':22s} {np.sqrt((r_yang**2).mean()):14.5f} "
      f"{np.sqrt((r_fit**2).mean()):14.5f}")
print(f"{'max |error|':22s} {np.abs(r_yang).max():14.5f} "
      f"{np.abs(r_fit).max():14.5f}")

# --- identifiability --------------------------------------------------------
# T_md(S) = T0 + cS*S is the density-maximum temperature.  In the polar box it
# lies far BELOW the freezing line, so the fit only ever samples one flank of
# the parabola and (Cb, T0, cS) are strongly degenerate -- only the local slope
# is constrained.  This must be stated, or the fitted values will be read as
# physically meaningful when they are not.
Tmd_lo = T0 + cS * SA_LO
Tmd_hi = T0 + cS * SA_HI
print()
print("-" * 74)
print("IDENTIFIABILITY")
print("-" * 74)
print(f"density-maximum temperature T_md(S) = T0 + cS*S:")
print(f"   at SA = {SA_LO}: {Tmd_lo:8.3f} C     at SA = {SA_HI}: {Tmd_hi:8.3f} C")
print(f"   freezing line spans           {gsw.CT_freezing(SA_HI, P_REF, 0.):.3f} "
      f"to {gsw.CT_freezing(SA_LO, P_REF, 0.):.3f} C")
if Tmd_hi < gsw.CT_freezing(SA_HI, P_REF, 0.0):
    print("   -> T_md is BELOW the freezing line across the whole box: the fit")
    print("      samples only ONE flank of the parabola, so Cb, T0 and cS are")
    print("      degenerate.  Only the local thermal expansion they imply is")
    print("      meaningful; do not read them as physical constants.")

try:
    J = fit.jac
    _, sv, _ = np.linalg.svd(J, full_matrices=False)
    print(f"   Jacobian condition number = {sv[0]/sv[-1]:.3e}"
          "   (large = degenerate, as expected)")
except Exception:
    pass

# --- what actually matters: the expansion coefficients ----------------------
print()
print("-" * 74)
print("THERMAL / HALINE EXPANSION AT THE BOX CENTRE  (what the flow feels)")
print("-" * 74)
Sc_, Tc_ = 0.5 * (SA_LO + SA_HI), -0.5
a_teos = gsw.alpha(Sc_, Tc_, P_REF)          # 1/K
b_teos = gsw.beta(Sc_, Tc_, P_REF)           # kg/g
rho_c = gsw.rho(Sc_, Tc_, P_REF)
# model: d rho/dT = -2*Cb*(T - T0 - cS*S) ; alpha = -(1/rho) d rho/dT
drhodT = -2 * Cb * (Tc_ - T0 - cS * Sc_)
drhodS = b0 + 2 * Cb * (Tc_ - T0 - cS * Sc_) * cS
a_fit = -drhodT / rho_c
b_fit = drhodS / rho_c
print(f"   at SA = {Sc_:.1f} g/kg, CT = {Tc_:.1f} C, p = {P_REF:.0f} dbar")
print(f"{'   alpha [1/K]':28s} TEOS-10 {a_teos:11.4e}   fit {a_fit:11.4e}"
      f"   ({100*(a_fit/a_teos-1):+.2f} %)")
print(f"{'   beta  [kg/g]':28s} TEOS-10 {b_teos:11.4e}   fit {b_fit:11.4e}"
      f"   ({100*(b_fit/b_teos-1):+.2f} %)")

# --- pressure sensitivity ---------------------------------------------------
print()
print("-" * 74)
print("PRESSURE SENSITIVITY (the simplified EOS carries no pressure term)")
print("-" * 74)
for p in (0.0, 50.0, 200.0):
    a = gsw.alpha(Sc_, Tc_, p)
    print(f"   p = {p:6.0f} dbar :  alpha = {a:.4e} 1/K "
          f"({100*(a/a_teos-1):+.2f} % vs surface)")
print("   -> thermobaricity is a real effect over a deep water column, but")
print("      negligible across the cm-to-tens-of-diameters DNS box, which is")
print("      the only vertical extent this model resolves.")

# --- nondimensional inputs --------------------------------------------------
print()
print("=" * 74)
print("NONDIMENSIONAL eos_* INPUTS")
print("=" * 74)
print("Convention (self-consistent, betaT = 1 by choice of the density scale):")
print("   theta = (T - T_i)/dT,   s = S/S_m,   drho_ref = Cb*dT^q")
print("   b = -|theta - Tmd0 - slope*s|^q + betaS*s")
print("   Tmd0  = (T0 - T_i)/dT")
print("   slope = cS*S_m/dT")
print("   betaS = b0*S_m/(Cb*dT^q)")
print()

CASES = [
    # name,            T_i,   T_ocean, S_m
    ("ice shelf / -1.9 C ice, 0.5 C ocean, S=34", -1.9, 0.5, 34.0),
    ("sea ice     / -1.9 C ice, 1.0 C ocean, S=32", -1.9, 1.0, 32.0),
]
q = 2.0
print(f"{'case':46s} {'dT':>6s} {'betaS':>10s} {'Tmd0':>9s} {'slope':>9s}")
for name, Ti, Tw, Sm in CASES:
    dT = Tw - Ti
    Tmd0 = (T0 - Ti) / dT
    slope = cS * Sm / dT
    betaS = b0 * Sm / (Cb * dT ** q)
    print(f"{name:46s} {dT:6.2f} {betaS:10.4f} {Tmd0:9.4f} {slope:9.4f}")

print()
print("Deck block for the first case:")
name, Ti, Tw, Sm = CASES[0]
dT = Tw - Ti
print("[eos]")
print(f"eos_q = {q}")
print("eos_betaT = 1.0")
print(f"eos_betaS = {b0 * Sm / (Cb * dT ** q):.6f}")
print(f"eos_Tmd0 = {(T0 - Ti) / dT:.6f}")
print(f"eos_Tmd_slope = {cS * Sm / dT:.6f}")
print()
print("VALIDITY RANGE (state this in the paper):")
print(f"   SA in [{SA_LO}, {SA_HI}] g/kg, CT in [{CT_LO}, {CT_HI}] C, "
      f"p ~ {P_REF} dbar,")
print(f"   RMS density error {np.sqrt((r_fit**2).mean()):.4f} kg/m^3, "
      f"max {np.abs(r_fit).max():.4f} kg/m^3.")
print("   No pressure/thermobaric terms.")

# --- the finding that actually matters for ECCO ------------------------------
print()
print("=" * 74)
print("SCALING CHECK -- the thermal density scale is the WRONG one for ECCO")
print("=" * 74)
print("betaS above is O(500-700).  That is not a bad fit; it is the physics:")
print("in polar water haline buoyancy dominates thermal buoyancy, the reverse")
print("of Yang's warm freshwater tank.  Quantifying it at the box centre:")
print()
for name, Ti, Tw, Sm in CASES:
    dT = Tw - Ti
    therm = a_teos * dT                       # fractional density change from dT
    for dS_label, dS in (("full meltwater contrast", Sm), ("1 g/kg stratification", 1.0)):
        hal = b_teos * dS
        print(f"   {name.split('/')[0].strip():10s} alpha*dT = {therm:.3e} vs "
              f"beta*dS = {hal:.3e}  ({dS_label}) -> haline/thermal = {hal/therm:7.1f}x")
print()
print("Consequence: with drho_ref = Cb*dT^q the buoyancy field is O(betaS) = O(700),")
print("so the 'free-fall' velocity built from that scale is wrong by ~sqrt(700) ~ 26x")
print("and Ra_T is not the governing number.  Nondimensionalize on the HALINE")
print("scale instead:  drho_ref = b0*S_m, which gives betaS = 1 and")
print("betaT = Cb*dT^q/(b0*S_m) = O(1e-3).  Same code, same equation, different")
print("input values -- the exponent q and the |.| form are untouched.")
print()
print(f"{'case':46s} {'betaT':>12s} {'betaS':>8s} {'Tmd0':>9s} {'slope':>9s}")
for name, Ti, Tw, Sm in CASES:
    dT = Tw - Ti
    print(f"{name:46s} {Cb*dT**q/(b0*Sm):12.6e} {1.0:8.1f} "
          f"{(T0-Ti)/dT:9.4f} {cS*Sm/dT:9.4f}")
print()
print("RECOMMENDED ECCO deck block (haline scaling, first case):")
name, Ti, Tw, Sm = CASES[0]
dT = Tw - Ti
print("[eos]")
print(f"eos_q = {q}")
print(f"eos_betaT = {Cb*dT**q/(b0*Sm):.6e}")
print("eos_betaS = 1.0")
print(f"eos_Tmd0 = {(T0-Ti)/dT:.6f}")
print(f"eos_Tmd_slope = {cS*Sm/dT:.6f}")
print()
print("With this scale Re = sqrt(Ra_S/Sc) and Ra_S = g*(b0*S_m/rho0)*H^3/(nu*kappa_S)")
print("is the governing Rayleigh number -- Ra_T becomes the subordinate one.")
