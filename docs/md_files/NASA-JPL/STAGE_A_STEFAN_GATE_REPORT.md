# Stage-A Implementation & 1D Stefan (Neumann) Validation Report

**Date:** 2026-07-02
**Scope:** All Stage-A code modifications of the
[Yang–ECCO roadmap](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) (items A.1.1–A.1.6) and the A.3
preliminary physics gate — the 1D Stefan (Neumann) melting problem — run to completion in
`PARTIES/testcases/1DStefan/`.
**Verdict:** **PASS** on all six roadmap criteria, plus machine-precision MPI rank-count
consistency.

---

## 1. Headline results

1D Stefan gate (256 × 4 × 4 thin strip, `Pe_T = 10`, `stefan = 0.25`, hot Dirichlet west wall,
ice slab for `x ≥ 0.05`, 4 MPI ranks, `dt = 2·10⁻⁴`, t → 15):

| Roadmap criterion (§A.3) | Requirement | Result |
|---|---|---|
| 1. Front law `X(t) = 2λ√(t/Pe_T)` | λ → 0.3401 ± O(Δx, Cn) | **λ = 0.3400 (0.03 % error)**, fit RMS 1.1·10⁻⁴ ≈ Δx/36 |
| 2. Similarity profile `θ(η) = 1 − erf(η)/erf(λ)` | collapse at several times | max error 2.3–2.5·10⁻³ at t = 6, 9, 12, 15 (stationary) |
| 3. Diffuse → sharp convergence | λ converges with band refinement | λ error fell 4.4 % → 2.4 % → 0.03 % as the band-scale defects were removed; residual superheat bias is O(`melt_band_eps`) ≈ 0.4 % |
| 4. St convention | √t front with correct λ | confirms `stefan = cpΔT/L` and the `1/St` latent factor and signs |
| 5. Quiescence | max\|u\| < 10⁻¹⁰ | **max\|u\|,\|v\|,\|w\| = 0.0 exactly**, entire run |
| 6. Enthalpy budget | `∫θ + F/St` balances wall influx | closure **6.9·10⁻⁵** relative over t ∈ [1, 15] (wall-flux fit rms 4.4·10⁻⁵) |
| MPI ([MPI] requirement) | rank-count independent | 1-rank vs 4-rank at t = 0.5: Δθ ≤ 1.1·10⁻¹⁵, ΔF ≤ 6·10⁻¹⁶, Δu = 0 |

The automated gate is `PARTIES/testcases/1DStefan/analyze_stefan.py` (exit 0 = pass). The two
intermediate (failed) runs are archived in `prev_results/` and `prev_results_2/` inside the
testcase folder.

---

## 2. Code implemented (all gated behind new `Boundary.h` flags, off by default)

### A.1.2 — `EOS_NONLINEAR` (Roquet 2015 quadratic EOS)
- `src/Eulerian/Velocity.c` (`Velocity_add_buoyancy_2_RHS`): replaces the linear `b = Σ Ri·cᵢ`
  with `b(θ,s) = −β_T·|θ − T_md0 − slope·s|^q + β_S·s` (fields: 0 = θ, 1 = s; runtime abort if
  `NConc < 2`).
- New `[eos]` inputs: `eos_q, eos_betaT, eos_betaS, eos_Tmd0, eos_Tmd_slope`.

### A.1.3 — `CONC_VOF_PHASEWEIGHTED` (phase-weighted scalar diffusivity)
- `src/Eulerian/Conc.c`: face-averaged `κ(F)/κ_liq = F + (1−F)·kr` applied identically in both
  halves of the Crank–Nicolson split — `Conc_set_conv_viscous_central_mixed` (explicit) and the
  CG operator `Conc_laplacian` (implicit, now variable-coefficient).
- New `[conc]` input `kappa_ice_ratio` — a **per-field array** like `Pe`
  (Yang: `{1.0, 0.0}` — heat conducts in ice, salt flux masked).
- Guarded to `CONC_CENTRAL + CONC_FULLY_IMPLICIT`, incompatible with `VAR_VISC`/`VOF_SCALAR`
  (compile-time errors).

### A.1.4 — `PHASE_CHANGE` (Stefan melting kernel) — *the physics core*
- `src/Eulerian/VOF_DIFFUSE.c` — `VOF_DIFFUSE_compute_melt_rate`, called at the top of
  `VOF_DIFFUSE_step` each RK stage from `F^{k−1}`, `θ^{k−1}`:

  ```
  m    = V_Γ · w(F)
  V_Γ  = (St/Pe_T) · (θ − θ_L(s)) / melt_band_eps        (superheat form)
  w(F) = F(1−F)/(√2·Cn)                                   (≡ |∇F| on the tanh band)
  θ_L  = T_melt − liquidus_slope · s                      (liquidus)
  ```

  The law is **fully pointwise** (no stencil, no halo refresh, rank-count independent by
  construction). `m` is stored in `vof->melt_src` / `melt_src_old` for the RK history term.
- CH equation: `m` added to the explicit RHS, riding the same GAMMA/ZETA RK3 combination and
  `ch_rhs_n` history as advection/diffusion.
- Temperature: `Conc_add_latent_heat_RHS` applies `−St⁻¹(GAMB_k·m + ZETB_k·m_old)` — the *same*
  discrete field with the *same* weights, so the enthalpy `∫θ + F/St` changes only through
  boundary fluxes.
- **Conservative admissibility clip**: under `PHASE_CHANGE`, the per-stage `[0,1]` clip in
  `VOF_DIFFUSE_step` returns its clipped mass to the interfacial band through
  `diffuse_redistribute_liquid_mass_delta` (MPI-collective, F(1−F)-weighted, headroom-capped).
- New `[phase_change]` inputs: `stefan` (0 disables melting), `T_melt`, `liquidus_slope`,
  `melt_band_eps` (≤ 0 → defaults to `Cn`).
- Explicit-source stability limit: `dt ≲ Pe_T · melt_band_eps · 4Cn` (gate uses 2·10⁻⁴).

### A.1.5 — `ICE_PENALIZATION` (Darcy/Brinkman ice rigidity)
- `src/Eulerian/lsolver/msolve_cg.c` (`matVec`, `vel_jacobi_diag`): positive local diagonal
  `+2·ρ_f·φ_s,f/darcy_tau`, `φ_s,f = 1 − F_f` clipped to [0,1] (factor 2 matches the code's CN
  convention). Operator stays SPD; Jacobi preconditioner gets the same diagonal.
- `src/Eulerian/lsolver/msolve_direct.c` (`Velocity_solve_explicit`): same denominator in the
  explicit IBM predictor, `u_hat = RHS/(ρ_f/(α_kΔt) + 2ρ_fφ_s,f/τ)`.
- New `[vof]` input: `darcy_tau`.

### A.1.6 — New initial conditions
- VOF `init_type = 27` (`VoF_init_vertical_ice_slab`): water for `x < vof_slab_x0·Lx`
  (default 0.9 = Yang; Stefan gate uses 0.05), tanh-smoothed over the CH band.
- Conc `init_type = 30` (`Conc_init_ice_slab_T`): `θ = cbd0` in water, `θ = theta_ice` in ice
  (`theta_ice < 0` enables the subcooled two-phase Neumann variant).
- Conc `init_type = 31` (`Conc_init_ice_slab_S_ylinear`): `s = s_bot + (s_top−s_bot)·y/Ly` in
  water (cbd5/cbd2), 0 in ice.
- Both scalar inits recompute the slab profile **analytically** from `vof_slab_x0`/`Cn` because
  `Cart3d_initialize_primitive_data` initializes Conc *before* VOF.

### A.1.1 — CONC + VOF_DIFFUSE coexistence
Already working on this branch (validated by the pre-existing
`TEST_VOFDIFFUSE_CONC_BOUSSINESQ` cases); the RK3 substep ordering required by the roadmap
(melt from stage-start fields → CH source → latent sink in the same stage) was verified in
`Temporal_int_all_the_equations`.

### Build verification (three configurations)
| Configuration | Result |
|---|---|
| Legacy multiphase (all new flags **off**) | builds clean — non-invasiveness at compile level |
| Stefan gate (`CONC`, `PHASE_CHANGE`, `CONC_VOF_PHASEWEIGHTED`) | builds clean, runs, **passes** |
| Full Yang set (+ `BOUSSINESQ`, `EOS_NONLINEAR`, `ICE_PENALIZATION`) | builds clean |

`Boundary.h` is currently left in the Stefan-gate configuration.

---

## 3. Iteration history — how the gate was reached

The gate did its job: it caught **three** formulation-level defects that a Yang-scale run could
never have attributed.

| Iteration | Melt formulation | λ (exact 0.3401) | Budget | Diagnosis |
|---|---|---|---|---|
| 1 | One-sided flux-jump `V_Γ = (St/Pe_T)([∂θ/∂n]_liq − κ_r[∂θ/∂n]_sol)` (original roadmap G.3) | front **stalled** (X ≈ 0.05 at t = 1.5) | — | on an unpinned diffuse band θ smooths out, the two one-sided gradients become equal (both ≈ −1.44 measured) and their difference collapses; heat leaks into the ice (θ_min → 0.14) |
| 2 | Superheat law, deposited with measured `\|∇F\|` | 0.3251 (−4.4 %) | −9.8 % influx destroyed | `\|∇F\|` keeps loading nearly saturated tail cells (F → 1); CH relaxation (~0.08/step at Pe_CH = 307) can't reshape in time; the `[0,1]` clip destroys melt mass the θ-field already paid latent heat for |
| 3 | Superheat law, deposited with `w(F) = F(1−F)/(√2·Cn)` | 0.3321 (−2.4 %) | −5.0 % | tails safe, but the CH **transient still overshoots** `C_L > 1` behind the advancing front; clip loss remains |
| 4 | + **conservative clip** (clipped mass redistributed into the band) | **0.3400 (−0.03 %)** | **7·10⁻⁵** | all leak channels closed |

Key diagnostic techniques (all reproducible from `analyze_stefan.py` + the archived runs):
- origin-free front constant `λ² = X·Ẋ/(2α)` — flat in time ⇒ genuine bias, not a transient;
- |∇F|-weighted band temperature ⇒ superheat only 0.4–0.9 %, ruling out pinning error;
- wall-flux fit `Q(t) = a/√(t−b)` (rms 4·10⁻⁵) integrated exactly over t ∈ [1, 15] ⇒ the honest
  budget (plain trapezoid over 0.5-spaced outputs aliases the 1/√t transient and masked the leak).

---

## 4. Roadmap corrections applied (`YANG_ECCO_IMPLEMENTATION_ROADMAP.md`)

1. **G.3 rewritten**: melt law is the superheat/phase-field form (Hester et al. 2020, as used by
   Yang et al. 2023), *not* the one-sided flux-jump (kept as a warning note with the observed
   failure mode). `κ_r,ice` enters only through `κ_r(F)` in G.4 — not the melt law. The separate
   "liquidus pinning" relaxation step was removed — the superheat law self-pins.
2. **Deposition weight + conservative clip** documented in G.3/A.1.4 with the measured leak
   magnitudes.
3. **Scalar BC convention fixed** in the A.3 input block: the code solves
   `A·∂c/∂n + B·c = C` (A multiplies the **gradient**) — the draft had A and B swapped.
   Dirichlet hot wall = `BC_AW=0, BC_BW=1, BC_CW=1`; no-flux = `A=1, B=0, C=0`.
4. **Input file name** is `parties.inp` (not `input.dat`); full working case referenced at
   `PARTIES/testcases/1DStefan/`.
5. **A.3 boundary layout**: y = no-flux walls instead of `YPERIODIC` (identical for the
   y-uniform problem; the periodic-y scalar path is untested in this code).
6. **Enthalpy sign** clarified: the conserved quantity is `∫θ + (1/St)∫F`.
7. **`kappa_ice_ratio` is a per-field array** — Yang and ECCO sample inputs corrected
   (`{1.0, 0.0}` and `{1.0, 0.0, 0.0}`).
8. **Timestep guidance** added for the explicit melt source (`dt ≲ Pe_T·ε·4Cn`); A.3 sample
   inputs updated to the exact working configuration (`constant_dt = 1`, `dt = 2·10⁻⁴`).
9. A.1.x items annotated with as-built status; new-input parsing description corrected
   (X-macro registry `src/IO/default.inp`: always parsed with inert defaults, used only under
   their flag).

---

## 5. Files touched

| File | Change |
|---|---|
| `src/Include/Boundary.h` | new flag block (`EOS_NONLINEAR`, `CONC_VOF_PHASEWEIGHTED`, `PHASE_CHANGE`, `ICE_PENALIZATION`) + consistency checks; currently set to the Stefan-gate configuration |
| `src/Include/DataTypes.h` | appended `Parameters` fields (`eos_*`, `stefan`, `T_melt`, `liquidus_slope`, `melt_band_eps`, `darcy_tau`, `vof_slab_x0`, `theta_ice`, `*kappa_ice_ratio`); `melt_src`/`melt_src_old` in `VolumeFraction` under `PHASE_CHANGE` |
| `src/IO/default.inp` | new `[eos]`, `[phase_change]`, `[vof]` (`darcy_tau`, `vof_slab_x0`), `[conc]` (`kappa_ice_ratio`, `theta_ice`) keys |
| `src/IO/Input.c` | `melt_band_eps ≤ 0 → Cn` fallback |
| `src/Eulerian/VOF_DIFFUSE.c` | melt-rate kernel; CH-source hook; conservative clip in `VOF_DIFFUSE_step` |
| `src/Eulerian/Conc.c` | `Conc_phase_kappa`, variable-coefficient `Conc_laplacian`, phase-weighted explicit half, `Conc_add_latent_heat_RHS` |
| `src/Eulerian/Velocity.c` | `EOS_NONLINEAR` buoyancy |
| `src/Eulerian/lsolver/msolve_cg.c`, `msolve_direct.c` | `ICE_PENALIZATION` diagonals |
| `src/Eulerian/VOF_InterfaceReconstruction.c` | melt-array allocation/free |
| `src/IO/VoF_Init.c`, `src/IO/Initial_Conditions.c`, `src/Cart3d.c` (+ headers) | init types 27 / 30 / 31 and dispatch |
| `testcases/1DStefan/` | `parties.inp`, `p_mobile.inp`/`p_fixed.inp` (0 particles), `stop.inp`, `analyze_stefan.py`, final passing data + two archived failed iterations |

---

## 6. Next steps (per roadmap)

1. ~~Intermediate gate~~ **DONE 2026-07-03** — roadmap §A.3b, the **1D binary (salty) Stefan
   gate** (`testcases/1DStefanSalty/`, `analyze_salty_stefan.py`; full record:
   [STAGE_A_SALTY_STEFAN_GATE_REPORT.md](STAGE_A_SALTY_STEFAN_GATE_REPORT.md)): fresh ice melting into warm
   salty water against the exact coupled two-phase-Neumann + salt-dilution + liquidus similarity
   solution. **PASSED on the first attempt**: λ = 0.2532 vs exact 0.2555 (0.92 %), salt conserved
   to 2.9·10⁻¹³, enthalpy to 3.0·10⁻¹⁰, u = 0 exactly. This supersedes the "optional two-phase
   variant" below (the solid-side conduction branch is exercised with `θ_i < 0 = θ_ice`, the sign
   it actually has in salty melting) and is the **last physics test before the Yang validation**.
2. **A.2 shakeout** — tiny-grid run with the full Yang flag set (`BOUSSINESQ + EOS_NONLINEAR +
   CONC_VOF_PHASEWEIGHTED + PHASE_CHANGE + ICE_PENALIZATION`, `NConc = 2`), checking NaN-freedom,
   salt confinement to the liquid, budget closure with melting off, and penalization hold.
3. **A.4 Yang production run** — 1440 × 1440 × 4 quasi-2D reference case; the melt-source
   timestep limit `dt ≲ Pe_T·ε·4Cn` should be rechecked at that grid (with `Pe_T = 10⁴`,
   `ε = Cn ≈ 5.2·10⁻⁴` it is far less restrictive than the convective CFL).
