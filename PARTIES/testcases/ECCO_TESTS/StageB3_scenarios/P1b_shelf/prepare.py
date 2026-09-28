#!/usr/bin/env python3
"""Generate every P1b input from one physical record; never submits jobs.

P1b (roadmap B.3, redesigned 2026-09-28): one fine-sand grain released beneath
an Antarctic ice-shelf base into the diffusive meltwater sublayer, melting ON.

Physical baseline (user decision 2026-09-28): 700 dbar, mixed layer below the
ice CT=0.0 C, SA=34.85 g/kg, friction velocity u*=0.5 cm/s; quartz grain
d=0.2 mm.  The interface state and sublayer thicknesses come from the
three-equation melt model (Jenkins et al. 2010 transfer coefficients).

Writes: case.json, parties.inp (production), p_mobile.inp, p_fixed.inp,
stop.inp, Boundary.scenario.h, gate/<variant>/parties.inp (1-D kernel gate),
gate/reference.json, smoke/parties.inp (full-size initial-condition check).
Needs numpy, scipy, gsw (Anvil: the anaconda 2025.06 python).
"""
import json
import math
from pathlib import Path
import re

import gsw
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erf, erfc

HERE = Path(__file__).resolve().parent
B2 = HERE.parents[1] / 'StageB2_closure' / 'B2_closure'
BASE = HERE.parent / 'P1a_feasibility' / 'parties.inp'

# ---------------------------------------------------------------- physics
P = 700.0                    # sea pressure, dbar (fixed-pressure EOS)
CT_M, SA_M = 0.0, 34.85      # mixed layer beneath the ice
USTAR = 0.005                # friction velocity, m/s
GAMMA_T, GAMMA_S = 0.011, 3.1e-4   # Jenkins, Nicholls & Corr (2010)
T_ICE_INTERIOR = -20.0       # only enters the three-equation heat balance
RHO_I, C_I = 917.0, 2009.0
D, RHO_P = 2.0e-4, 2650.0    # fine quartz sand
NU, KT, D_PHYS, SC = 1.8e-6, 1.4e-7, 1.0e-9, 70.0
CP, LATENT = 3991.86795711963, 333500.0   # TEOS-10 cp0 for CT; constant L
G = 9.81

# ---------------------------------------------------------------- numerics
N_PER_D = 24
LX, LY, LZ = 8.0, 22.0, 8.0
Y_INT = 20.0                 # ice occupies y in [20, 22]
GAP = 0.35                   # grain top below y_int, clear of the tanh band
T_END, DT_FIELD, DT_MAX, DT_START = 30.0, 0.1, 0.005, 0.0005


def rho(sa, ct):
    return gsw.rho_t_exact(sa, gsw.t_from_CT(sa, ct, P), P)


def metrics(e):
    e = np.asarray(e)
    return dict(max_absolute=float(np.max(np.abs(e))), rms=float(np.sqrt(np.mean(e * e))))


def three_equation():
    rho_w = float(rho(SA_M, CT_M))

    def resid(sb):
        tb = gsw.CT_freezing(sb, P, 0)
        wb = USTAR * GAMMA_S * (SA_M - sb) / (RHO_I / rho_w * sb)
        return (rho_w * 3974.0 * USTAR * GAMMA_T * (CT_M - tb)
                - RHO_I * wb * (LATENT + C_I * (tb - T_ICE_INTERIOR)))
    sb = brentq(resid, 1.0, SA_M - 1e-9)
    tb = float(gsw.CT_freezing(sb, P, 0))
    wb = USTAR * GAMMA_S * (SA_M - sb) / (RHO_I / rho_w * sb)
    return float(sb), tb, float(wb)


def eos_fit(t_ref, dt, s_ref):
    """Existing EOS_NONLINEAR q=2 form, fitted as in P1a: density-maximum
    locus pinned to GSW, curvature from d(rho)/dCT, then haline slope and a
    constant (absorbed by pressure)."""
    s_lo, s_hi, t_lo, t_hi = 10.0, 36.0, -2.2, 0.4
    sa, ct = [a.ravel() for a in np.meshgrid(np.linspace(s_lo, s_hi, 131),
                                             np.linspace(t_lo, t_hi, 105))]
    dens = rho(sa, ct)
    smd = np.linspace(s_lo, s_hi, 53)
    tmd = np.array([minimize_scalar(lambda t: -float(rho(s, t)), bounds=(-12., 6.),
                                    method='bounded').x for s in smd])
    k, t0 = np.polyfit(smd, tmd, 1)
    eps = 1e-4
    drdt = (rho(sa, ct + eps) - rho(sa, ct - eps)) / (2 * eps)
    x = -2 * (ct - t0 - k * sa)
    b = float(np.dot(x, drdt) / np.dot(x, x))
    a, h = np.linalg.lstsq(np.array([np.ones(sa.size), sa]).T,
                           dens + b * (ct - t0 - k * sa) ** 2, rcond=None)[0]
    assert b > 0 and h > 0
    sv, tv = np.meshgrid(np.linspace(s_lo, s_hi, 197), np.linspace(t_lo, t_hi, 157))
    exact = rho(sv, tv)
    fitted = a - b * (tv - t0 - k * sv) ** 2 + h * sv
    dt_exact = (rho(sv, tv + eps) - rho(sv, tv - eps)) / (2 * eps)
    thermal_err = -2 * b * (tv - t0 - k * sv) - dt_exact
    liquid = tv >= gsw.CT_freezing(sv, P, 0)
    coeffs = dict(eos_q=2.0, eos_betaT=float(b * dt ** 2 / (h * s_ref)), eos_betaS=1.0,
                  eos_Tmd0=float((t0 - t_ref) / dt), eos_Tmd_slope=float(k * s_ref / dt))
    record = dict(choice='EOS_NONLINEAR q=2, density-maximum locus pinned to GSW at 700 dbar',
                  gsw_version=gsw.__version__, pressure_dbar=P,
                  density_routine='rho_t_exact(SA, t_from_CT(SA,CT,p), p)',
                  SA_range_g_kg=[s_lo, s_hi], CT_range_C=[t_lo, t_hi],
                  dimensional_a_b_Tmd0_TmdSlope_h=[float(a), b, float(t0), float(k), float(h)],
                  density_error_kg_m3=metrics((fitted - exact)[liquid]),
                  density_error_over_haline_scale=metrics((fitted - exact)[liquid] / (h * s_ref)),
                  thermal_derivative_error_kg_m3_K=metrics(thermal_err[liquid]),
                  density_maximum_locus_error_K=metrics(t0 + k * smd - tmd),
                  note='Errors over liquid (above-freezing) states of the fit box. The interface water '
                       'sits near its density maximum at this pressure, which a linear EOS cannot represent.',
                  coefficients=coeffs)
    return coeffs, record, float(h)


def stefan_similarity(pe_t, pe_s, st, tm, lam_liq):
    """Binary Stefan similarity: semi-infinite ice at theta=0, water theta=1,
    s=1, equal diffusivity in both phases, salt only in the liquid."""
    le = pe_s / pe_t

    def g(x):
        return math.sqrt(math.pi) * x * math.exp(x * x) * (1 + math.erf(x))

    def resid(lam):
        si = 1 / (1 + g(lam * math.sqrt(le)))
        thi = tm - lam_liq * si
        return (math.sqrt(math.pi) * lam * math.exp(lam * lam) / st
                - ((1 - thi) / (1 + math.erf(lam)) - thi / math.erfc(lam)))
    lam = brentq(resid, 1e-8, 2.0)
    si = 1 / (1 + g(lam * math.sqrt(le)))
    return dict(lambda_=lam, s_interface=si, theta_interface=tm - lam_liq * si,
                front=f'X(t) = 2*lambda*sqrt(t/Pe_T), Pe_T={pe_t}')


def write_deck(path, base, overrides, extra_conc=None, header=''):
    text = '\n'.join(l for l in base.splitlines() if not l.startswith('#')) + '\n'
    for key, value in overrides.items():
        text, n = re.subn(r'^' + re.escape(key) + r'\s*=.*$', f'{key} = {value}', text, flags=re.M)
        assert n == 1, (key, n)
    if extra_conc:
        block = '\n'.join(f'{k} = {v}' for k, v in extra_conc.items())
        text = text.replace('[conc]\n', '[conc]\n' + block + '\n', 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(header + text)


def main():
    sb, tb, wb = three_equation()
    t_ref, dt = tb, CT_M - tb                 # theta = (T - T_b)/(T_M - T_b)
    s_ref = SA_M                               # s = SA/SA_M
    delta_T, delta_S = KT / (USTAR * GAMMA_T), D_PHYS / (USTAR * GAMMA_S)
    ell_T, ell_S = 2 * delta_T / math.sqrt(math.pi) / D, 2 * delta_S / math.sqrt(math.pi) / D
    s_i = sb / SA_M

    eos, eos_record, h = eos_fit(t_ref, dt, s_ref)
    rho0 = float(rho(SA_M, CT_M))
    beta = h / rho0
    uref = math.sqrt(G * beta * s_ref * D)
    tref, re_, pr = D / uref, uref * D / NU, NU / KT
    pe_t, pe_s, st = re_ * pr, re_ * SC, CP * dt / LATENT
    gstar, rho_s = 1 / (beta * s_ref), RHO_P / rho0
    ga = re_ * math.sqrt((rho_s - 1) * gstar)
    rep = brentq(lambda r: r * (1 + 0.15 * r ** 0.687) - ga ** 2 / 18, 1e-4, 1e3)

    # Liquidus: secant through the interface state and the far water, so the
    # initial interface is exactly on the liquidus (theta_L(s_i) = 0).
    tf_m = float(gsw.CT_freezing(SA_M, P, 0))
    m_f = (tb - tf_m) / (SA_M - sb)
    a_f = tb + m_f * sb
    sgrid = np.linspace(10, 35, 251)
    liq_err = a_f - m_f * sgrid - gsw.CT_freezing(sgrid, P, 0)
    t_melt, lam_liq = (a_f - t_ref) / dt, m_f * s_ref / dt
    assert abs(t_melt - lam_liq * s_i) < 1e-12

    cn = 0.75 / N_PER_D
    relax = 4 * pe_t * cn * 2 * cn
    ch_dt_limit = 2.5 / (0.5 * 12 * N_PER_D ** 2 / (0.9 / cn))
    nx, ny, nz = int(LX * N_PER_D), int(LY * N_PER_D), int(LZ * N_PER_D)
    cells = nx * ny * nz
    y_grain = Y_INT - GAP - 0.5
    su_step_p1a = 182.7555556 / 5032 / 2211840         # measured, per cell-step
    steps = T_END / DT_MAX
    su = steps * cells * su_step_p1a

    # initial profile table (dimensional), grain-relevant numbers
    zs = np.array([0, 0.5, 1, 2, 3, 4, 6, 8, 10, 15, 20])
    T_prof = tb + dt * erf(zs / ell_T)
    S_prof = sb + (SA_M - sb) * erf(zs / ell_S)
    rho_prof = rho(S_prof, T_prof) - rho0
    z = np.linspace(0, 20, 4001)
    rz = rho(sb + (SA_M - sb) * erf(z / ell_S), tb + dt * erf(z / ell_T))
    n2 = G / rho0 * np.gradient(rz, z * D)
    nmax = float(np.sqrt(n2.max()))
    ws = rep * NU / D

    gate = dict()
    for name, n in [('G1_sliq1_n24', 24), ('G0_sliq0_n24', 24),
                    ('G1_sliq1_n24_halfdt', 24), ('G1_sliq1_n32', 32)]:
        gate[name] = n
    ref = stefan_similarity(pe_t, pe_s, st, t_melt, lam_liq)

    record = dict(
        design_date='2026-09-28',
        decision='User: Antarctic shelf, fine sand, free grain in the meltwater sublayer, melting ON; '
                 'baseline CT_M=0 C and u*=0.5 cm/s (2026-09-28).',
        physical=dict(p_dbar=P, CT_mixed_layer_C=CT_M, SA_mixed_layer_g_kg=SA_M, ustar_m_s=USTAR,
                      Gamma_T=GAMMA_T, Gamma_S=GAMMA_S,
                      Gamma_source='Jenkins, Nicholls & Corr (2010) JPO 40, 2298',
                      T_ice_interior_C=T_ICE_INTERIOR, d_m=D, rho_p_kg_m3=RHO_P,
                      nu_m2_s=NU, kappa_T_m2_s=KT, salt_diffusivity_physical_m2_s=D_PHYS,
                      Sc_sim=SC, cp_J_kg_K=CP, latent_J_kg=LATENT, rho0_kg_m3=rho0,
                      note='Idealized constants. The three-equation state sets the interface and sublayer '
                           'thicknesses; the box itself is quiescent (slack-current limit).'),
        three_equation=dict(SA_interface_g_kg=sb, CT_interface_C=tb, melt_rate_m_s=wb,
                            melt_rate_m_yr=wb * 3.156e7, delta_T_m=delta_T, delta_S_m=delta_S,
                            delta_T_d=delta_T / D, delta_S_d=delta_S / D,
                            meltwater_fraction_interface=1 - sb / SA_M,
                            caveat='S_b depends on Gamma_T/Gamma_S; Keitzl et al. (2016) argue the '
                                   'three-equation form overestimates interface salinity by up to 40%.'),
        eos=eos_record,
        liquidus=dict(a_f_C=a_f, m_f_K_per_g_kg=m_f, T_melt=t_melt, liquidus_slope=lam_liq,
                      error_over_SA_10_35_K=metrics(liq_err),
                      note='Secant through (S_b, T_b) and (S_M, T_f(S_M)) at 700 dbar; theta_L(s_i)=0 exactly.'),
        scaling=dict(T_ref_C=t_ref, DeltaT_K=dt, S_ref_g_kg=s_ref, beta_per_g_kg=beta,
                     U_ref_m_s=uref, t_ref_s=tref, Re=re_, Pr=pr, Pe_T=pe_t, Pe_S=pe_s, Le=pe_s / pe_t,
                     St=st, G_star=gstar, rho_s=rho_s, Ga=ga, Re_p_SchillerNaumann=rep,
                     w_s_m_s=ws, w_s_over_Uref=ws / uref, N_max_s=nmax, Froude_ws_over_Nd=ws / (nmax * D)),
        initial_state=dict(theta='erf(zeta/ell_T)', s='F*[s_i+(1-s_i)*erf(zeta/ell_S)]',
                           tracer='F*(1-s_i)*(1-erf(zeta/ell_S))', ell_T_d=ell_T, ell_S_d=ell_S,
                           s_interface=s_i, ell_definition='ell = 2*delta/sqrt(pi): wall gradient equals '
                           'the three-equation flux',
                           profile_z_d=zs.tolist(), CT_C=T_prof.tolist(), SA_g_kg=S_prof.tolist(),
                           rho_minus_rho_M_kg_m3=rho_prof.tolist()),
        geometry=dict(L_d=[LX, LY, LZ], N=[nx, ny, nz], cells=cells, cells_per_d=N_PER_D,
                      y_interface=Y_INT, ice_thickness_d=LY - Y_INT, grain_center=[LX / 2, y_grain, LZ / 2],
                      grain_top_gap_d=GAP, lateral_period_debris_volume_fraction=math.pi / 6 / LX ** 2,
                      bottom='no-slip wall, insulated; top: ice, insulated'),
        numerics=dict(Cn=cn, melt_band_eps=2 * cn, Pe_CH=0.9 / cn, darcy_tau=1e-3,
                      max_dt=DT_MAX, default_dt=DT_START, cfl=0.3, T_END=T_END, T_END_s=T_END * tref,
                      field_interval=DT_FIELD, melt_relaxation_time=relax,
                      ch_explicit_dt_limit=ch_dt_limit, liquid_referenced_salinity=1,
                      yang_salt_transport=0, F_release=-1.0),
        cost=dict(su_per_cell_step_P1a=su_step_p1a, nominal_steps=steps, nominal_SU=su,
                  SU_with_30pct=1.3 * su, storage_GB_per_frame=cells * (526102296 + 235839912) / 2211840 / 1e9,
                  frames=int(T_END / DT_FIELD) + 1),
        gate=dict(reference=ref, variants=list(gate),
                  design='Quiescent 1-D strip (4 x 20d x 4 cells-wide), ice on top at theta=0 (s=0), water '
                         'theta=1, s=1, EOS buoyancy off; same Pe_T, Pe_S, St, liquidus, band as P1b.'),
        references=['https://journals.ametsoc.org/view/journals/phoc/40/10/2010jpo4317.1.xml',
                    'https://arxiv.org/abs/1606.03004',
                    'https://teos-10.github.io/GSW-Python/density.html'])
    (HERE / 'case.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')

    base = BASE.read_text()
    phys = dict(Re=re_, stefan=st, T_melt=t_melt, liquidus_slope=lam_liq,
                melt_band_eps=2 * cn, Cn=cn, Pe_CH=0.9 / cn,
                meltwater_tracer=1, yang_salt_transport=0, liquid_referenced_salinity=1,
                Pe='{' + ', '.join(map(repr, [pe_t, pe_s, pe_s])) + '}',
                rho_s=rho_s, rho_prImp=rho_s, grav='{0.0, ' + repr(-gstar) + ', 0.0}',
                F_release=-1.0, darcy_tau=1e-3, cfl=0.3,
                BC_AN='{1.0, 1.0, 1.0}', BC_BN='{0.0, 0.0, 0.0}', BC_CN='{0.0, 0.0, 0.0}',
                BC_AS='{1.0, 1.0, 1.0}', BC_BS='{0.0, 0.0, 0.0}', BC_CS='{0.0, 0.0, 0.0}',
                **eos)
    prod = dict(xmax=LX, ymax=LY, zmax=LZ, NXM=nx, NYM=ny, NZM=nz,
                time_max=T_END, output_time_interval=DT_FIELD, output_time_interval_2d=DT_FIELD,
                default_dt=DT_START, max_dt=DT_MAX, vof_slab_x0=repr(Y_INT / LY),
                conc_init_type='{36, 37, 38}', theta_ice=0.0, cbd0=1.0, cbd5=1.0, cbd2=1.0, **phys)
    sub = dict(sublayer_ell_T=repr(ell_T), sublayer_ell_S=repr(ell_S), sublayer_s_interface=repr(s_i))
    write_deck(HERE / 'parties.inp', base, prod, sub,
               '# P1b: released fine-sand grain beneath an Antarctic ice-shelf base. See case.json.\n')
    (HERE / 'p_mobile.inp').write_text(f'1\n{LX/2} {y_grain} {LZ/2} 0.5\n')
    (HERE / 'p_fixed.inp').write_text('0\n')
    (HERE / 'stop.inp').write_text('0\n')
    (HERE / 'Boundary.scenario.h').write_bytes((B2 / 'Boundary.validation.h').read_bytes())

    # Matched ice-only control (no grain): same deck, coarser field cadence.
    # Difference fields P1b - control isolate what the grain does to the layer.
    ctrl = dict(prod, output_time_interval=0.5, output_time_interval_2d=0.5)
    write_deck(HERE / 'ctrl/parties.inp', base, ctrl, sub,
               '# P1b ice-only control (no grain). See case.json.\n')
    (HERE / 'ctrl/p_mobile.inp').write_text('0\n')

    smoke = dict(prod, time_max=0.2, output_time_interval=0.1, output_time_interval_2d=0.1)
    write_deck(HERE / 'smoke/parties.inp', base, smoke, sub,
               '# P1b smoke: full-size deck, initial condition + first steps only.\n')

    # 1-D kernel gate: thin strip, sharp step initial condition, buoyancy off.
    for name, n in gate.items():
        cng = 0.75 / n
        w = 4.0 / n
        g = dict(phys, xmax=repr(w), ymax=20.0, zmax=repr(w), NXM=4, NYM=int(20 * n), NZM=4,
                 time_max=80.0, output_time_interval=2.0, output_time_interval_2d=80.0,
                 default_dt=DT_START / (2 if 'halfdt' in name else 1),
                 max_dt=DT_MAX / (2 if 'halfdt' in name else 1),
                 vof_slab_x0=0.5, Cn=repr(cng), melt_band_eps=repr(2 * cng), Pe_CH=repr(0.9 / cng),
                 conc_init_type='{34, 35, 0}', theta_ice=0.0, cbd0=1.0, cbd2=1.0, cbd5=1.0,
                 eos_betaT=0.0, eos_betaS=0.0,
                 liquid_referenced_salinity=0 if 'sliq0' in name else 1)
        write_deck(HERE / 'gate' / name / 'parties.inp', base, g, None,
                   f'# P1b 1-D kernel gate {name}. See case.json gate block.\n')
    (HERE / 'gate' / 'reference.json').write_text(json.dumps(dict(ref, Pe_T=pe_t, Pe_S=pe_s, St=st,
                                                                  T_melt=t_melt, liquidus_slope=lam_liq,
                                                                  y_interface=10.0), indent=2) + '\n')
    print(json.dumps(dict(three_equation=record['three_equation'], scaling=record['scaling'],
                          liquidus=record['liquidus'], eos_errors=dict(
                              density=eos_record['density_error_over_haline_scale'],
                              thermal=eos_record['thermal_derivative_error_kg_m3_K']),
                          coefficients=eos, gate_reference=ref, cost=record['cost'],
                          initial_state={k: record['initial_state'][k] for k in ['ell_T_d', 'ell_S_d', 's_interface']},
                          numerics=record['numerics']), indent=1))


if __name__ == '__main__':
    main()
