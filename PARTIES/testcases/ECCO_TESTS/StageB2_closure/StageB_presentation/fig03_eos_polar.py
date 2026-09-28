#!/usr/bin/env python3
"""Figure 3 — re-fitting the equation of state to polar conditions (B.1.1 / B.1.1a).

Panels (a)-(c) are recomputed here with gsw, using the same functional form,
the same polar box and the same two scenarios as
``ECCO_TESTS/StageB1_development/eos_polar_fit.py``; running this script reproduces that script's
printed table.  Panel (d) plots the recorded B.1.1a extrapolation table.

The headline is not the fit quality.  It is the scaling consequence: in polar
water haline buoyancy dominates thermal buoyancy by two orders of magnitude, so
an ECCO deck must nondimensionalize on the HALINE density scale.  Carrying
Yang's thermal scale over gives betaS ~ 691 and a free-fall velocity wrong by
sqrt(691) ~ 26x.
"""
import json
import sys
from pathlib import Path

import gsw
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stageb_common as C

plt = C.use_style()
CM = C.cmaps()
N = C.ROADMAP_NUMBERS

# --- the polar box, and the fit (same construction as eos_polar_fit.py) ----
SA_LO, SA_HI, CT_LO, CT_HI, P_REF = 28.0, 35.0, -2.2, 2.0, 0.0
NG = 60
SAg, CTg = np.meshgrid(np.linspace(SA_LO, SA_HI, NG),
                       np.linspace(CT_LO, CT_HI, NG), indexing='ij')
mask = CTg >= gsw.CT_freezing(SAg, P_REF, 0.0)
s_pt, t_pt = SAg[mask], CTg[mask]
rho_ref = gsw.rho(s_pt, t_pt, P_REF)


def rho_model(par, S, T, q=2.0):
    rho0, Cb, T0, cS, b0 = par
    return rho0 - Cb * np.abs(T - T0 - cS * S) ** q + b0 * S


YANG = np.array([1000.0, 0.011, 4.0, -0.25, 0.77])
fit = least_squares(lambda p: rho_model(p, s_pt, t_pt) - rho_ref, YANG,
                    method='lm', max_nfev=200000)
rho0, Cb, T0, cS, b0 = fit.x
err_fit = fit.fun
err_yang = rho_model(YANG, s_pt, t_pt) - rho_ref
rms_fit, rms_yang = np.sqrt((err_fit ** 2).mean()), np.sqrt((err_yang ** 2).mean())

# scenarios, exactly as in eos_polar_fit.py
CASES = [('ice shelf', -1.9, 0.5, 34.0), ('sea ice', -1.9, 1.0, 32.0)]
# alpha and beta are evaluated once at the BOX CENTRE, as eos_polar_fit.py does,
# so the two cases differ only through their dT and S_m.
S_C, T_C = 0.5 * (SA_LO + SA_HI), -0.5
alpha = float(gsw.alpha(S_C, T_C, P_REF))
beta = float(gsw.beta(S_C, T_C, P_REF))
rows = []
for name, T_i, T_o, S_m in CASES:
    dT = T_o - T_i
    betaS_thermal = b0 * S_m / (Cb * dT ** 2)
    betaT_haline = Cb * dT ** 2 / (b0 * S_m)
    rows.append(dict(name=name, dT=dT, S_m=S_m, betaS_thermal=betaS_thermal,
                     betaT_haline=betaT_haline,
                     aT=alpha * dT, bS_full=beta * S_m, bS_1=beta * 1.0,
                     Tmd0=(T0 - T_i) / dT, slope=cS * S_m / dT))

(C.DATADIR / 'eos_polar_fit.json').write_text(json.dumps(dict(
    box=dict(SA=[SA_LO, SA_HI], CT=[CT_LO, CT_HI], p=P_REF),
    yang=dict(zip(['rho0', 'Cb', 'T0', 'cS', 'b0'], YANG.tolist()), rms=rms_yang),
    polar=dict(zip(['rho0', 'Cb', 'T0', 'cS', 'b0'], fit.x.tolist()), rms=rms_fit),
    cases=rows), indent=2) + '\n')

fig = plt.figure(figsize=(12.6, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.26,
                      left=0.068, right=0.975, top=0.845, bottom=0.095)

# --- (a) fit error over the polar box ------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 0]))
ax.grid(False)
E = np.full(SAg.shape, np.nan)
E[mask] = np.abs(err_fit)
pm = ax.pcolormesh(SAg, CTg, E, cmap=CM['blue'], shading='auto',
                   vmin=0, vmax=np.nanmax(E))
cb = fig.colorbar(pm, ax=ax, pad=0.02)
cb.set_label('|error| vs TEOS-10  (kg m$^{-3}$)', fontsize=8.5, color=C.INK2)
cb.ax.tick_params(labelsize=8, colors=C.MUTED)
cb.outline.set_edgecolor(C.AXIS)
sline = np.linspace(SA_LO, SA_HI, 200)
ax.plot(sline, gsw.CT_freezing(sline, P_REF, 0.0), color=C.STATUS['critical'],
        lw=1.6, label='freezing line')
ax.set_xlabel('absolute salinity  $S_A$  (g kg$^{-1}$)')
ax.set_ylabel('conservative temperature  $\\Theta$  (°C)')
ax.set_title('a   The re-fit is 770× better inside the polar box')
ax.legend(loc='upper right')
ax.text(SA_LO + 0.12, CT_HI - 0.36,
        f'RMS error\n  Yang lab coefficients   {rms_yang:.4f}\n'
        f'  polar re-fit                   {rms_fit:.4f}',
        fontsize=8.5, color=C.INK2, va='top', linespacing=1.6)

# --- (b) the crossover ---------------------------------------------------
ax = C.finish(fig.add_subplot(gs[0, 1]))
S = np.linspace(0.05, 36, 400)
Tf = gsw.CT_freezing(S, P_REF, 0.0)


def t_maxdensity(SA, p=0.0):
    T = np.linspace(-15.0, 12.0, 4001)
    r = gsw.rho(np.full_like(T, SA), T, p)
    return float(T[int(np.argmax(r))])


Tmd = np.array([t_maxdensity(v) for v in S])
cross = float(np.interp(0.0, (Tmd - Tf)[::-1], S[::-1]))
ax.fill_between(S, Tf, 12, color='#eef3fa', lw=0)
ax.plot(S, Tf, color=C.BLUE, lw=2.0, label='freezing point  $T_f(S)$')
ax.plot(S, Tmd, color=C.ORANGE, lw=2.0, label='density maximum  $T_{md}(S)$')
ax.axvline(cross, color=C.MUTED, lw=1.2, ls='--')
ax.text(cross - 0.5, 10.4, f'crossover  S ≈ {cross:.0f}', color=C.INK2,
        fontsize=8.5, ha='left', rotation=90, va='top')
ax.axvspan(SA_LO, 36, color='#fdf1e9', lw=0, zorder=0)
ax.text(31.5, 4.5, 'ECCO', color=C.ORANGE, fontsize=10, ha='center',
        fontweight='semibold')
ax.axvline(5.0, color=C.AQUA, lw=1.6)
ax.text(5.6, 4.5, 'Yang, $S_m$ = 5', color=C.AQUA, fontsize=9)
ax.text(34, 8.6, 'liquid', color=C.MUTED, fontsize=9, ha='right')
ax.set_xlim(0, 36)
ax.set_ylim(-6.6, 11)
ax.set_xlabel('salinity  $S$  (g kg$^{-1}$)')
ax.set_ylabel('temperature  (°C)')
ax.set_title('b   Stage A validated the machinery, not the regime')
ax.legend(loc='upper left')
ax.text(0.6, -5.7,
        'beyond the crossover the density maximum is below the freezing\n'
        'line, so meltwater is unambiguously buoyant',
        fontsize=8.2, color=C.INK2, linespacing=1.5, va='bottom')

# --- (c) haline vs thermal buoyancy --------------------------------------
# A dot plot, not bars: on a log axis a bar has no meaningful baseline, and
# bars drawn from the axis floor would also bury the legend.
ax = C.finish(fig.add_subplot(gs[1, 0]))
ax.set_xscale('log')
ypos = np.arange(len(rows))[::-1]
SERIES = [('thermal   $\\alpha\\,\\Delta T$', 'aT', C.ORANGE, 'o'),
          ('haline   $\\beta\\,\\Delta S$   full meltwater',
           'bS_full', C.BLUE, 'o'),
          ('haline   $\\beta\\,\\Delta S$   1 g kg$^{-1}$',
           'bS_1', '#86b6ef', 's')]
for yv, r in zip(ypos, rows):
    ax.plot([r['aT'], r['bS_full']], [yv, yv], color=C.AXIS, lw=2.5,
            solid_capstyle='round', zorder=2)
for label, key, col, mk in SERIES:
    ax.plot([r[key] for r in rows], ypos, mk, ms=10, color=col, mec=C.SURFACE,
            mew=1.4, ls='none', zorder=4, label=label)
for yv, r in zip(ypos, rows):
    ax.text(np.sqrt(r['aT'] * r['bS_full']), yv + 0.17,
            f"{r['bS_full']/r['aT']:.0f}×", ha='center', fontsize=11,
            color=C.STATUS['critical'], fontweight='semibold')
    ax.text(r['bS_1'], yv - 0.22, f"{r['bS_1']/r['aT']:.1f}×", ha='center',
            va='top', fontsize=9, color=C.INK2)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{r['name']}\nΔT = {r['dT']} K,  S = {r['S_m']:.0f} g kg$^{{-1}}$"
                    for r in rows], fontsize=9)
ax.set_xlim(3e-5, 3e-1)
ax.set_ylim(-1.35, ypos[0] + 0.75)
ax.set_xlabel('buoyancy contrast (dimensionless)')
ax.set_title('c   In polar water salt, not heat, drives the buoyancy')
ax.grid(axis='y', visible=False)
ax.legend(loc='lower left', fontsize=8.4)
ax.text(2.9e-1, -1.28,
        'So an ECCO deck must scale on $\\Delta\\rho = b_0 S_m$:\n'
        f"betaS = 1, betaT = {rows[0]['betaT_haline']:.3e} (ice shelf).\n"
        f"Yang's thermal scale gives betaS = {rows[0]['betaS_thermal']:.0f} —\n"
        'a free-fall velocity wrong by √691 ≈ 26×.',
        fontsize=8.2, color=C.STATUS['critical'], va='bottom', ha='right',
        linespacing=1.5)

# --- (d) extrapolation, from the recorded B.1.1a table -------------------
ax = C.finish(fig.add_subplot(gs[1, 1]))
d = N['dilution']
Sv = np.array(d['S'], float)
ax.axhline(0, color=C.AXIS, lw=1.0)
ax.axvspan(SA_LO, 36, color='#f0efec', lw=0, zorder=0)
ax.text(35.6, 0.02, 'fit box', color=C.MUTED, fontsize=8, ha='right')
ax.plot(Sv, d['quadratic'], 'o-', color=C.ORANGE, ms=6,
        label='Roquet quadratic, vertex pinned')
ax.plot(Sv, d['linear'], 's-', color=C.BLUE, ms=6, label='linear')
ax.annotate(f"{d['quadratic'][-1]:+.3f}", (Sv[-1], d['quadratic'][-1]),
            textcoords='offset points', xytext=(10, -2), fontsize=9,
            color=C.ORANGE, va='center')
ax.annotate(f"{d['linear'][-1]:+.3f}", (Sv[-1], d['linear'][-1]),
            textcoords='offset points', xytext=(10, 2), fontsize=9,
            color=C.BLUE, va='center')
ax.set_xlim(-2, 36)
ax.invert_xaxis()
ax.set_xlabel('salinity the meltwater plume visits  (g kg$^{-1}$)   →  fresher')
ax.set_ylabel('density-anomaly error vs TEOS-10  (kg m$^{-3}$)')
ax.set_title('d   Outside the box the linear model wins by ≈ 5×')
ax.set_ylim(-0.40, 0.055)
ax.legend(loc='lower left')
ax.text(34.5, -0.10,
        'inside the box both are ~100× more accurate than needed\n'
        f"(rms/std: linear {N['eos_inbox']['linear  c0 + c1θ + c2s']} %, "
        f"quadratic {N['eos_inbox']['Roquet quadratic, vertex pinned']} %)\n"
        'recommendation: for ECCO do not compile EOS_NONLINEAR',
        fontsize=8.2, color=C.INK2, linespacing=1.5, va='top')

fig.suptitle('Stage B · the equation of state, re-fitted for polar water',
             x=0.068, ha='left', fontsize=15, fontweight='semibold', color=C.INK, y=0.965)
fig.text(0.068, 0.895,
         'B.1.1 / B.1.1a — input values only, no code change. The fitted Cb, T0, cS are '
         'individually degenerate (Jacobian condition number 3.5e4); quote α and β, not them.',
         ha='left', fontsize=9.5, color=C.INK2)

C.save(fig, 'fig03_eos_polar',
       'Panels a–c recomputed with gsw against TEOS-10, same construction as '
       'ECCO_TESTS/StageB1_development/eos_polar_fit.py (fit written to data/eos_polar_fit.json). '
       'Panel d is the recorded B.1.1a extrapolation table. ΔT and S_m are the two '
       'placeholder scenarios; the nondimensional block must be regenerated once the '
       'ECCO case is fixed.')
