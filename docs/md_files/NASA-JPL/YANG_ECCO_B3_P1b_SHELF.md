# B.3 P1b — released fine-sand grain beneath an Antarctic ice-shelf base

**Status (2026-10-01):** **production segment complete and analysed (§8.7,
[RESULTS_20977443.md](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf/RESULTS_20977443.md)).**

**Status (2026-09-30, 21:50):** all three production jobs are running and healthy. The grain is at
terminal settling and carries about one grain volume of meltwater down (§8.6).

**Status (2026-09-30, 16:50):** all three production jobs are running and healthy (§8.5).
The earlier status follows.

**Status (2026-09-30, 11:30):**
- The first P1b (**20973340**, Yang operator) **failed at t = 0.738**: the Yang salt field
  blew up at the moving grain (§8.4).
- P1b was resubmitted with the **default operator** as **20977443**, with a matched control
  **20977444**.
- The Yang control **20973341** keeps running as an operator-comparison reference.

This file is the controlling
record for P1b. The roadmap
[§B.3](YANG_ECCO_IMPLEMENTATION_ROADMAP.md#b3-full-ecco-production-inputs--scenario-design-and-two-precursors)
summarizes it and links here.

Inputs, generator and analysis scripts:
[`StageB3_scenarios/P1b_shelf/`](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf/).
Output: `/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/`.

---

## 1. Why P1 was redesigned

The first P1 ([P1a, job 20906327](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/RESULTS_20906327.md))
waited for a locked 1 mm grain to melt out under 0 °C water at p = 0. It showed
slight net freezing, no release and a failed relative rigidity check. The
2026-09-28 investigation found three separate problems.

1. **The freezing was a diffuse-interface artifact, not the physics of those
   inputs.** The exact binary-Stefan solution for P1a's inputs *melts*
   (λ = 0.0112, interface θ = 0.046, s = 0.954). A 1-D enthalpy model of the finite
   column, cold top included, keeps melting (0.37 d lost by t = 2×10⁴). P1a ran with
   `liquid_referenced_salinity = 0`: the liquidus then reads the volume-averaged
   salinity `F·s`, about half the liquid value at mid-band. The resulting temperature
   error is roughly 0.5 × `liquidus_slope`. That is 0.007 for Yang and B.2
   (slope 0.014) but 0.5 for P1a (slope 1.01). The band pinned at θ ≈ 0.44 instead
   of 0.05, above the θ = 0.2 threshold beyond which the thin cold ice removes more
   heat than the water supplies. The `F∇s` operator also drove `s/F` up to 3.2 on the
   ice side of the band. The negative tracer came from the resulting spurious freezing.
2. **Waiting for melt-out is unaffordable for any realistic polar forcing.** A front
   rise of 0.709 d is needed for the shell trigger. At d = 1 mm this takes
   3–6×10⁴ code units (0.5–1 h physical), about 0.5–1×10⁶ SU at P1a's timestep.
   The explicit Cahn–Hilliard term caps the timestep near 0.021 on this grid
   (Yang production ran at 96% of the same limit), which gives only a 10× saving.
   The grain then falls through the box in about 0.8 code units.
3. **Release timing is not a deliverable anyway.** B.2 showed escape timing is
   ill-conditioned and does not converge in dt.

The rigidity "failure" was a ratio of two negligible speeds in a still box
(domain maximum 2 µm/s).

**Consequence.** Release becomes the **initial condition**, as B.3.1 already
allowed (`F_release = -1`). P1b studies what happens after the grain leaves the ice:
- settling through the meltwater sublayer;
- meltwater carried down by the grain;
- wake mixing;
- the interplay with continued melting (melting stays ON, by user decision).

## 2. Physical case (user decisions 2026-09-28)

**Antarctic ice-shelf base, fine sand, melting ON.**

| Quantity | Value | Source / note |
|---|---|---|
| Sea pressure | 700 dbar | warm-cavity shelf base; fixed-pressure EOS |
| Mixed layer beneath the ice | CT = 0.0 °C, SA = 34.85 g/kg | deep water diluted by meltwater, 2.4 K above freezing |
| Friction velocity u* | 0.5 cm/s | currents ≈ 0.1 m/s; sets sublayer thickness and melt rate |
| Transfer coefficients | Γ_T = 0.011, Γ_S = 3.1×10⁻⁴ | Jenkins, Nicholls & Corr (2010) |
| Grain | quartz, d = 0.2 mm, 2650 kg/m³ | fine sand, the dominant released fraction |
| Ice | fresh glacial ice, isothermal at the interface temperature, flat | basal gradient ≈ 0.06 K/mm, so 0.02 K across the block |

**Three-equation interface state:**

| Quantity | Value |
|---|---|
| Interface salinity / temperature | S_b = 21.30 g/kg, T_b = −1.670 °C |
| Meltwater fraction at the ice | 39% |
| Melt rate | 35.0 m/yr (Pine Island main shelf observes 10–30 m/yr) |
| δ_T = κ/(u*Γ_T) | 2.55 mm = 12.7 d |
| δ_S = D/(u*Γ_S) | 0.645 mm = 3.2 d |

The interface salinity depends only on the water temperature and the Γ ratio; u*
sets the thicknesses and the melt rate. Keitzl et al. (2016) argue the
three-equation form overestimates interface salinity by up to 40%, so a fresher
variant is a natural follow-up.

**Initial profile below the ice** (z = distance below the interface; ρ relative
to the mixed layer at 700 dbar):

| z | CT (°C) | SA (g/kg) | meltwater fraction | ρ − ρ_M (kg/m³) |
|---|---|---|---|---|
| 0 | −1.670 | 21.30 | 0.389 | −10.750 |
| 0.5 d | −1.605 | 23.39 | 0.329 | −9.079 |
| 1 d | −1.539 | 25.40 | 0.271 | −7.471 |
| 2 d | −1.409 | 28.93 | 0.170 | −4.649 |
| 4 d | −1.158 | 33.22 | 0.047 | −1.225 |
| 6 d | −0.926 | 34.58 | 0.008 | −0.150 |
| 10 d | −0.542 | 34.85 | 0 | +0.037 |
| 20 d | −0.082 | 34.85 | 0 | +0.006 |

Salt controls the density: the layer is light down to about 7 d. Below it sits a
thin cold, salty zone slightly heavier than ambient (+0.04 kg/m³ near 9 d),
because the heat sublayer is 4× thicker than the salt one. This is a real
double-diffusive feature; it evolves over about 40 s, against a 0.82 s window.
The interface water lies near its density maximum at 700 dbar, so the thermal
expansion coefficient varies about 100× across the layer. A linear EOS cannot
represent that.

**Grain physics:**

| Quantity | Value |
|---|---|
| Ga / Re_p (Schiller–Naumann) | 6.17 / 1.73 |
| Settling speed w_s | 15.6 mm/s |
| Froude number w_s/(N d) at the steepest gradient | 6.2 |
| Layer density jump vs grain excess density | 0.66% |
| Crossing the light layer | ≈ 3 code units |

Expect a modest effect of the layer on the grain. The main signal is the reverse:
meltwater carried down by the grain, and its wake.

## 3. Model setup

**Scaling** (B.3.2 haline form, length = d):

| Quantity | Value |
|---|---|
| θ | (T − T_b)/(T_M − T_b), ΔT = 1.670 K |
| s | SA/34.85 |
| U_ref / t_ref | 7.30 mm/s / 27.4 ms |
| Re / Pe_T / Pe_S (Sc = 70) | 0.812 / 10.43 / 56.81 |
| St | 0.01999 |
| G* | 36.78 |
| ρ_s | 2.5699 |

**EOS:** `EOS_NONLINEAR` q = 2, refitted at 700 dbar over SA 10–36, CT −2.2…+0.4:
- density-maximum locus pinned to GSW, curvature fitted to dρ/dCT, then haline
  slope and constant;
- liquid-state errors: density ≤ 0.044 kg/m³ (0.16% of the haline scale),
  dρ/dCT ≤ 0.0021 kg m⁻³ K⁻¹;
- coefficients: `eos_betaT = 7.0319e-4`, `eos_Tmd0 = 2.54472`,
  `eos_Tmd_slope = −4.75943`, `eos_betaS = 1`;
- `richardson = {0,0,0}`.

**Liquidus:** secant through (S_b, T_b) and (S_M, T_f(S_M)) at 700 dbar:
- `T_melt = 0.727676`, `liquidus_slope = 1.190446`;
- θ_L(s_b) = 0 exactly;
- error ≤ 0.028 K over 10–35 g/kg.

**Diffuse band — melting ON:**
- **Yang finite-interface salt operator** (`yang_salt_transport = 1`,
  `liquid_referenced_salinity = 0`), decided 2026-09-30 after two gate rounds (§8.2).
- The salt field is the liquid salinity S. It evolves with the flux (F+δ)∇S plus the
  meltwater-dilution term −S·ṁ/(F+δ), and its conserved quantity is ∫(F+δ)S.
- Init 37 writes S unmasked, also inside the ice, so the band starts in equilibrium:
  S = s_i at θ = 0 on the liquidus.
- The tracer keeps the volume-averaged convention.

**Numerics:**
- band: Cn = 0.03125, `melt_band_eps` = 0.0625, Pe_CH = 28.8, `darcy_tau` = 1e-3
  with the hard threshold (B.2-validated);
- timestep: `max_dt` = 0.005, CFL 0.3, `default_dt` = 5e-4. Limits:

  | Limit | Value |
  |---|---|
  | Cahn–Hilliard explicit | 0.021 |
  | Melt relaxation time | 0.082 |
  | CFL at w_s | ≈ 0.0058 |

**Ice block:**
- VOF init 28, interface at y = 20, ice from y = 20 to 22 (2 d);
- θ = 0 (interface temperature) throughout;
- salinity and tracer zero inside the ice;
- top boundary insulated for all scalars.

**Meltwater layer** (new Conc init types, ζ = max(y_int − y, 0)):

```
T (36):      θ = erf(ζ/ℓ_T)                              ℓ_T = 2δ_T/√π = 14.361 d
S (37):      s = F·[s_i + (1 − s_i)·erf(ζ/ℓ_S)]           ℓ_S = 2δ_S/√π = 3.640 d, s_i = 0.611264
tracer (38): c = F·(1 − s_i)·[1 − erf(ζ/ℓ_S)]             meltwater fraction w.r.t. the mixed layer
```

- ℓ = 2δ/√π makes the wall gradient equal to the three-equation flux.
- Salt and tracer share a diffusivity, so c = 1 − s holds in the liquid for all time,
  and the melt source adds meltwater with c = 1, s = 0.
- δ_S uses the physical salt diffusivity; Sc = 70 broadens the layer by only about
  5% over the window.

**Grain:**
- centre (4, 19.15, 4), radius 0.5, starting at rest;
- top 0.35 d below the ice interface, clear of the tanh band (1 − F = 4×10⁻⁴ there);
- `F_release = -1`;
- starts in water that is about 30% meltwater.

**Box:**
- 8 × 22 × 8 d at 24 cells/d, giving 192 × 528 × 192 = 19.46 M cells;
- periodic in x and z;
- no-slip, insulated bottom;
- the lateral period corresponds to one grain per 8 × 8 d, about 0.8% by volume for a
  one-grain-thick debris sheet (dispersed basal debris ice is 1–10% by mass).

**Duration and output:**
- T_END = 30 (0.82 s physical): the grain lands around t ≈ 9, and the displaced
  layer then restratifies (buoyancy period ≈ 18 code units);
- fields every 0.1 (301 frames, about 6.7 GB each, about 2 TB total); profiles every
  10 steps.

**Ice-only control:** the same deck without the grain, fields every 0.5.
P1b − control isolates the grain's effect (roadmap B.3.5: attribution needs a
matched background).

## 4. Code change (additive, no solver change)

`src/IO/Initial_Conditions.c` gains `Conc_init_ice_sublayer()`. `Cart3d.c`
dispatches Conc init types **36/37/38** to it.

New `[conc]` keys `sublayer_ell_T`, `sublayer_ell_S` and `sublayer_s_interface`
are registered in `default.inp` with inert defaults, and appended to `Parameters`.
The new init types reject a non-positive ℓ or a non-positive `cbd5`. Existing init
types and decks are untouched, per the B.1.0 contract (additive init cases,
appended fields).

The build is isolated at `P1b_shelf/build/`, made by `build_validation.py
--configuration sediment` with the header identical to
`B2_closure/Boundary.validation.h`.

## 5. Pre-flight: 1-D kernel gate and smoke test

**Gate.** It tests the band's melt/liquidus/salt kernel in P1b's regime
(liquidus slope 1.19, Le = 5.44, St = 0.02), which no earlier gate reached.

Setup:
- a quiescent strip, 4 × 20 d × 4 cells;
- ice above y = 10 at θ = 0, s = 0; water θ = 1, s = 1;
- buoyancy off, same band and timestep as P1b.

Exact binary-Stefan reference: **λ = 0.019382**, interface s = 0.92215 and
θ = −0.37010.

Four variants:
- `G1_sliq1_n24`: the P1b setting;
- `G0_sliq0_n24`: P1a's default band;
- `G1_sliq1_n24_halfdt`;
- `G1_sliq1_n32`: P2 resolution.

Criteria (declared in `analyze_gate.py` before any result):

| Criterion | Requirement |
|---|---|
| C1 | \|λ_fit/λ − 1\| ≤ 5% |
| C2 | no net freezing at any profile time |
| C3 | interface θ and s/F within 0.05 of exact |
| C4 | enthalpy ≤ 10⁻³ of the latent heat spent; salt ≤ 10⁻⁸; tracer − melt ≤ 10⁻³ |
| C5 | n32 within 3% and half-dt within 1% of n24 |

**P1b is submitted only if a band variant passes C1–C4.**

**Smoke test.** The full production deck to t = 0.2 checks:
- the initial profiles against the formulas;
- ice and grain geometry;
- first-step melt direction;
- grain acceleration toward w_s;
- measured cost per step.

## 6. What P1b measures, and its acceptance

**Measures:**
- grain velocity against height through and below the sublayer;
- meltwater tracer carried below the initial layer, relative to the control;
- wake tracer filament and its decay;
- layer restratification and dissipation;
- the melt rate with a falling grain present.

**Health criteria:**

| Check | Requirement |
|---|---|
| Salt drift | ∫(F+δ)S drift reported. The 1-D gate shows 1.8×10⁻⁵ over t = 80 with no grain; the Yang form is not exactly conservative (§8.2) |
| Tracer vs net melt | ≤ 10⁻³ once melt is significant |
| Heat budget | B.2 normalization ≤ 10⁻⁴ (insulated box) |
| Fields | finite |
| Liquid tracer | ≥ −10⁻⁶ |
| Net ice change | not negative, and no new ice below the band |
| Bulk-ice speed | ≤ 1% of the domain maximum; meaningful now that the grain moves |

**Interpretation limits:**
- the box is still water (slack current);
- the sublayer is a snapshot and is not maintained by turbulence;
- periodic neighbours hinder settling more at Re ≈ 2 than in B.2;
- no terminal-velocity or ocean flux-law claim;
- attribution of mixing to the grain uses the control.

## 7. Cost and next steps

| Run | Cells | Nominal steps | SU (P1a throughput) | With 30% |
|---|---|---|---|---|
| P1b | 19.46 M | 6,000 | 1,918 | 2,493 |
| Ice-only control | 19.46 M | 6,000 | ≈ 1,900 | ≈ 2,500 |
| P2 (32/d, later) | 46.1 M | ≈ 8,000 | ≈ 6,000 | ≈ 7,800 |

Allocation before submission: 71,516.5 SU. P2 waits for the P1b/control analysis.

## 8. Run record

### 8.1 Pre-flight round 1 — debug job 20950025 (2026-09-28): gate FAILS, smoke FAILS

About 44 SU (256 cores × 10.4 min). Output:
`P1_Attempts/P1b_shelf/gate_smoke_20950025/`; verdict in `gate/gate_result.json`.

**Gate:**

| Variant | λ_fit (exact 0.019382) | Error | Interface θ (exact −0.370) | Interface s/F (exact 0.922) | max s/F, ice side |
|---|---|---|---|---|---|
| `G1_sliq1_n24` (P1b setting) | 0.03758 | **+93.9%** | −1.051 | 1.262 | 2.7 |
| `G1_sliq1_n24_halfdt` | 0.03758 | +93.9% | −1.051 | 1.262 | 2.7 |
| `G1_sliq1_n32` | 0.03558 | +83.6% | −0.932 | 1.261 | 2.3 |
| `G0_sliq0_n24` (P1a default) | 0.01760 | **−9.2%** | −0.224 | 1.641 | **140** |

- **What passes:** all variants conserve exactly (enthalpy ≤ 1.5×10⁻⁷ of the latent heat
  spent, salt ≤ 1.5×10⁻¹⁴, tracer − melt ≤ 3×10⁻⁵), and none shows net freezing
  (C2, C4 pass).
- **What fails:** C1 and C3 fail everywhere.
- **It is not a numerical-resolution problem:** halving dt changes nothing, and 32/d moves
  λ by only −5%. It is the model form.
- **Mechanism:** the default salt operator ∇·(F∇s) acts on volume-averaged salt, whose
  equilibrium is s = constant. It therefore drives salt into the partly frozen band. The
  interface's liquid salinity reaches 1.26–1.64× the far field, where the exact interface
  is diluted to 0.92.
  - With liquid referencing, the enrichment lowers the liquidus on the ice side and doubles
    the melt rate.
  - Without it, the band reads volume-averaged salinity: λ is 9% slow with the wrong
    interface state.
  - At Yang's liquidus slope (0.014) both errors were invisible. At 1.19 they are not.
- **Decision:** per the pre-declared rule, P1b is not submitted with either default-operator
  variant. The code's Yang finite-interface operator (`yang_salt_transport = 1`, flux
  (F+δ)∇S on the liquid salinity, plus the meltwater dilution term) is the physically
  consistent alternative. It becomes gate round 2, with the same criteria written into
  `analyze_gate.py` before any round-2 result.
- **Supporting code change (additive):** init 37 writes the liquid salinity unmasked when
  the Yang operator is on.

**Smoke** (full size, 224 ranks):
- the first RK stage needed two pressure passes (divergence 1.0×10⁻³ → 4×10⁻⁷);
- in the second stage, the grain's first motion, the divergence stalled at **2.82×10⁻⁴**,
  so the run aborted "Poisson solver did not converge" (loop cap 10 passes, tolerance
  10⁻⁶);
- B.2's moving-grain run always converged in one pass (≤ 2.4×10⁻⁷ over 25,215 steps);
- the uneven rank split is not the cause: `NX = NXM+1` makes every PARTIES split uneven,
  P1a included;
- the projection is a standard constant-density update, and the boundary routine does not
  overwrite interior velocity.

Seven one-factor diagnostics on a small 4 × 13 × 4 d box localize it, using a diagnostic
binary that prints the divergence location (`build_diag`, never used for production):
grain distance from the ice, no grain, no buoyancy, Re × 100, Yang operator, melting off.

### 8.2 Pre-flight round 2 — debug job 20971700 (2026-09-29): cause found, Yang operator adopted

Run on 1 node (128 ranks) for 8.2 min, about 17.5 SU. It replaced the never-started 2-node
job 20970943: in the debug queue, 1-node jobs start while 2-node requests wait on
resources. Output: `P1_Attempts/P1b_shelf/preflight2_20971700/`.

**Poisson stall — root cause found.** Diagnostics used a 4 × 13 × 4 d box, the first ~17
steps, and a build that prints the divergence location.

| Run | Change | Result | Divergence after projection |
|---|---|---|---|
| D1 | none (liquid-referenced band) | **stalls** | 8.56×10⁻⁵ at (42, 233, 42), inside the grain |
| D2 | grain 3 d below the ice | **stalls** | 1.82×10⁻⁴ at (53, 177, 37), inside the grain |
| D5 | Re × 100 | **stalls** | — |
| D7 | melting off | **stalls** | — |
| D3 | no grain | runs | — |
| D4 | buoyancy off | runs | ~10⁻¹¹ |
| D6 | Yang operator (liquid referencing off) | runs | ~10⁻¹¹; salt CG ≤ 18 iterations; grain accelerates normally |

**Cause.** `liquid_referenced_salinity = 1` makes the EOS divide the salt by the liquid
fraction C_L. Inside a resolved grain C_L → 0 (the sediment occupies the cell).
- At t = 0, s/(C_L+δ) reaches **1,564** in the grain cells; physical values are ≤ 1.
- The immersed boundary forces those cells while a buoyancy about 10³× too large acts
  on them.
- The pressure solve then works on enormous magnitudes. The stalled divergences are
  quantized at 2⁻²⁷, which is round-off from velocities of order 10⁶.
- It is **not** the uneven rank split, low Re, the ice, or melting.

**Rule:** do not use `liquid_referenced_salinity = 1` together with resolved grains.

**Gate round 2 (Yang operator).** Reference λ = 0.019382, interface s = 0.92215 and
θ = −0.37010.

| Variant | λ_fit | Error | θ error | s error | Enthalpy / latent | Salt ∫F·S drift | Tracer − melt |
|---|---|---|---|---|---|---|---|
| `G2_yang_n24` | 0.019413 | **+0.16%** | +0.017 | +0.001 | 2.4×10⁻⁷ | **1.8×10⁻⁵** | 1.4×10⁻⁵ |
| `G2_yang_n32` | 0.019405 | +0.12% | +0.013 | +0.0008 | 8.7×10⁻⁸ | 1.5×10⁻⁵ | 1.6×10⁻⁵ |

C1, C2, C3 and C5 pass. The ice side of the band stays at S ≤ 1, so there is no salt
pile-up. **C4 fails on salt only:** 1.8×10⁻⁵ against the pre-declared 10⁻⁸. Strictly,
`analyze_gate.py` therefore still reports "STOP".

**Decision (2026-09-30, recorded as a deviation from the pre-declared rule).** P1b uses
the Yang operator.
- **Why the bound was mis-specified:** it came from the default operator, which conserves
  ∫s exactly. The Yang form evolves the liquid salinity non-conservatively. ∫(F+δ)S is
  preserved only as far as the liquid fraction changes through the melt term; the
  Cahn–Hilliard relaxation and the clip/restore also reshape it. Stage A measured the
  same order for this operator.
- **Size:** 1.8×10⁻⁵ of total salt is about 0.05% of P1b's meltwater salt deficit
  (about 4% of the box's salt).
- **Accuracy:** it is the only variant that reproduces the exact solution: 0.2% on λ,
  against −9% and +94% for the default operator.
- **The only alternative band setting is already ruled out:** liquid referencing breaks
  the pressure solve with a resolved grain.
- **Still unmeasured:** the salt budget with a **moving** grain under this operator is new.
  P1b's own profiles measure ∫F·S every 10 steps.

### 8.3 Production submission — 2026-09-30

| Job | Case | Resources | Status |
|---|---|---|---|
| **20973340** | P1b (grain) | 1 node, 128 ranks, `shared`, 24 h | **running** since 2026-09-30 01:06 |
| **20973341** | ice-only control | same | **running** since 2026-09-30 01:06 |

- Binary: `build_v2`, SHA-256 `27c2cdcf…d761f`.
- Staged from `P1_Attempts/P1b_shelf/submit_P1b_20260930/`; runs write to
  `run_prod_<job>` and `run_ctrl_<job>`.
- Balance before submission: 71,454.5 SU.
- Hard cap: 2 × 3,072 SU. Nominal: about 1.9k SU each.
- Checks to run on the first frames: `check_smoke.py` against Data_0, the Poisson loop
  converging in one pass, and the ∫F·S drift from the profiles.

**Early health check (P1b at t = 0.13–0.23; control at t = 0.19–0.31; both at 20–30 min):**
- **Initial state:** `check_smoke.py` on P1b's Data_0 matches θ, S (unmasked, Yang), the
  tracer and the ice tanh to ≤ 2.2×10⁻¹⁶. The tracer minimum is 0.
- **Grain:** its support spans y = 18.69–19.65. It is moving, y 19.150 → 19.081 by
  t = 0.129, at v = −0.845 (0.40 w_s) and still accelerating. That is plausible at
  Re ≈ 2, 0.35 d from the ice ceiling, where wall drag slows it.
- **Pressure:** always **one pass**. Maximum divergence after projection is 8.6×10⁻¹¹
  (P1b) and 3.7×10⁻¹³ (control); there are no NaNs or errors.
- **Budgets:**

  | Budget | P1b | Control |
  |---|---|---|
  | Tracer − melt | ~10⁻⁹ | ~10⁻⁹ |
  | Liquid salt ∫F·S drift | −8.8×10⁻⁷ after 0.07 d of grain travel | +7.8×10⁻¹⁰ |

  The P1b drift is the new moving-grain term (the profiles exclude C_S ≥ 0.05 cells,
  whose membership moves with the grain). Its growth over the transit will be reported,
  not extrapolated.
- **Melting:** −1.5×10⁻³ d³ of ice in 0.1 code units, about 41 m/yr. The initial profile
  was built to give the three-equation 35 m/yr.

**Throughput is below plan.** P1b takes about 24.7 s per step (that sample included one
frame write), and the control about 20 s, against the ~9 s implied by P1a's throughput.

| Solve (per call) | P1b | P1a |
|---|---|---|
| HYPRE pressure | 5.7 s, 57 iterations | 0.17 s, 34 iterations |
| Velocity CG | 14–16 iterations | ≈ 4 iterations |
| Yang salt | 8.5 iterations | 1.3 iterations (default operator) |

The pressure solve dominates. It has 8.8× more cells per rank, and Re ≈ 0.8 stiffens the
viscous solves.

**Consequences for this segment:**
- the 24 h wall limit's checkpoint-stop (23.5 h) should leave P1b near t ≈ 18 and the
  control near t ≈ 21;
- each costs its 3,072 SU hard cap;
- the grain should land around t ≈ 9–11, so the segment covers the transit and about
  8 code units of wake/restratification;
- finishing to T_END = 30 would cost roughly 2.0k SU (P1b) plus 1.3k SU (control) by
  restart;
- that is **not automatic**: it will be decided after the segment is analyzed.

### 8.4 P1b-Yang failure (job 20973340) and the switch to the default operator

**What happened.**
- 20973340 started 01:06 and failed at 02:15 (1 h 09 min, 147 SU) at step 178, t = 0.738,
  with "Poisson solver did not converge". The residual stuck at 1.07×10⁻⁶ against the 10⁻⁶
  tolerance.
- The grain had reached v = −1.49 (0.70 w_s) at y = 18.33. It then decelerated, and in the
  last step v jumped from −1.39 to −0.98 while its lateral velocity grew 20–40×.
- The saved frames show why: **the Yang salt field diverged**.

  | t | S in fluid cells | Buoyancy minimum | Ice below y = 19 |
  |---|---|---|---|
  | 0.2 | [0.611, 1.000] | — | 1×10⁻¹¹ |
  | 0.5 | [−2.75, 4.00] | −2.9 | 1.4×10⁻⁷ |
  | 0.6 | [−27, 59] | −39 | 7×10⁻⁵ |
  | 0.7 | [−427, 237] | −3,331 | 5×10⁻³ |

  The spurious ice comes from the wild salinity shifting the local liquidus.
- The first out-of-range cells, at t = 0.3–0.4 (S = 0.49, then 0.13), sit on the grain's
  **trailing (top) surface**: sediment-support cells with C_S = 0.05–0.10 and C_L = 0.

**Cause.**
- The Yang operator scales each cell's diffusion by 1/(F+δ) and uses face capacity
  δ + ½(F₁+F₂).
- A support cell (F = 0) beside water (F ≈ 0.95) therefore gets a coefficient of about
  0.47/δ ≈ 5×10⁵.
- The Crank–Nicolson split does not damp such stiff modes, and the moving grain changes
  those coefficients every step.
- The 1-D gate (a smooth tanh band, no jump) and the no-grain control contain no such jump,
  so both stay stable. The Yang control **20973341** was healthy at t = 6.9 after 10 h.
- The volume-averaged **tracer**, which uses the *default* operator in the same run, stayed
  bounded ([−6×10⁻²⁰, 0.35]) throughout. That is direct evidence the default operator is
  stable with this moving grain.

**Decision (2026-09-30).**
- P1b and its matched control run with the **default operator, no liquid referencing**:
  - it is B.2-validated for moving grains and conserves ∫s exactly;
  - the known cost is band-local (gate round 1, G0): λ −9% and a wrong interface salinity;
  - the t = 0 band is not in equilibrium under the mixture liquidus. Expect a brief freezing
    transient that pins the interface near θ ≈ 0.36 (about 0.6 K too warm), with small
    negative tracer confined to the band.
  - these artifacts stay within about 0.26 d of the ice, which the grain leaves by t ≈ 0.2.
- The Yang control (20973341) continues. Comparing it with the default control (20977444)
  **measures directly how much the band-kernel choice changes the meltwater sublayer** over
  the window. That decides whether P1b-default needs a Yang rerun later.
- **Fixing the Yang operator at resolved sediment** is a development item (for example,
  capacity max(F, C_S)+δ in support cells). It needs a moving-grain test before use.

| Job | Case | Operator | Status |
|---|---|---|---|
| 20973340 | P1b | Yang | **failed** t = 0.738, 147 SU |
| 20973341 | control | Yang | running (reference) |
| **20977443** | P1b | default | submitted 2026-09-30 |
| **20977444** | control | default | submitted 2026-09-30 |

### 8.5 Mid-run health check — 2026-09-30, 16:50

| Job | Case | Runtime | Reached | Rate | Projected at the 23.5 h stop |
|---|---|---|---|---|---|
| **20977443** | P1b, default operator | 4.2 h | t = 3.18 | 0.76 t.u./h | t ≈ 17.7 |
| **20977444** | control, default operator | 5.0 h | t = 4.22 | 0.84 t.u./h | t ≈ 19.8 |
| 20973341 | control, Yang operator | 15.8 h | t = 10.5 | 0.67 t.u./h | t ≈ 15.6 |

No run has had a pressure-loop failure.

**P1b-default is numerically healthy past the Yang failure point (t = 0.738).**
- **Grain:** falling steadily with lateral drift of ~10⁻⁶. At t = 3.18 it is at y = 13.80,
  6.2 d below the interface, with v = −1.941.

  | t | v / w_s |
  |---|---|
  | 0.5 | 0.64 |
  | 1.0 | 0.78 |
  | 2.0 | 0.88 |
  | 3.2 | 0.91 |

  Extrapolating the fall suggests a landing around t ≈ 10.
- **Budgets:**
  - salt ∫s drift ≤ 3.3×10⁻¹⁵ (exact);
  - tracer − melt ≤ 1.9×10⁻⁷;
  - heat drift 1.7×10⁻⁶, against 18 units of latent heat.
- **Health:**
  - liquid (F > 0.95) tracer stays ≥ 0;
  - liquid salinity stays ≤ 1;
  - bulk-ice speed is 0.41% (t = 1), 0.10% and 0.02% of the domain maximum, inside the
    1% criterion.

**Expected default-kernel artifacts, now measured.**
1. **Net freezing, independent of the grain.** P1b and the default control both lost
   −0.353 d³ by t = 3; the rate is slowing (−0.19 at t = 0.65, −0.28 at 1.5). The Yang
   control **melts**, +0.063 d³. The three-equation rate would give about +0.03 d³. The
   volume-averaged liquidus starts the band undercooled, so the **default kernel gets the
   sign of melting wrong here**, not just the rate.
2. **Wake ice, P1b only.** 13,196 cells below y = 19.5 hold more than 10⁻³ ice, up to 2.4%,
   at y = 18.5–19.5 and within 0.84 d of the grain's path (median 0.44). The default control
   has none. Probable mechanism: the wake pulls the band's diffuse tail down, and the
   undercooled volume-averaged liquidus then freezes it (compare the B.2 "ice tongue" note).
3. **Sublayer contamination is confined near the ice.** Comparing the two controls (no
   grain) at t = 3, default minus Yang:

   | z below ice | liquid S diff | θ diff |
   |---|---|---|
   | 0.10 d | −0.094 | +0.199 |
   | 0.19 d | −0.104 | +0.179 |
   | 0.48 d | −0.032 | +0.118 |
   | 0.73 d | −0.007 | +0.075 |
   | 0.98 d | −0.001 | +0.044 |
   | 1.48 d | 0.000 | +0.011 |
   | ≥ 1.98 d | 0.000 | ≤ 0.002 |

   The tracer profiles agree everywhere. The thermal difference spreads as √t: θ is
   0.24 (0.4 K) too warm at the interface.

**Interpretation for P1b's science.**
- Below about 2 d from the ice the default and Yang states are identical at t = 3; the grain
  was already 2.1 d below the interface at t = 1. The settling through the lower sublayer
  and the mixed layer, the wake below about 2 d, and the grain-carried tracer are therefore
  not affected by the kernel choice at this level.
- Flagged as model-contaminated:
  - the top ~1.5–2 d, growing as √t;
  - the net melt sign;
  - the wake ice near the ice.
- Melt-rate claims need the Yang operator repaired at resolved sediment.

**Cost outlook.**
- All three production runs will reach their 3,072-SU caps, about 9.2k SU.
- The campaign total is about 9.4k SU with the failed Yang P1b (147) and the pre-flights (62).
- That leaves roughly 62k of the 71.5k SU available at the start.

### 8.6 Progress check — 2026-09-30, 21:50

| Job | Case | Reached | Runtime | Rate | Notes |
|---|---|---|---|---|---|
| 20977443 | P1b, default operator | t = 7.26 | 9.2 h | 0.79 t.u./h | 73 frames, 451 GB; projected stop at t ≈ 18.6, ~12:10 Oct 1 |
| 20977444 | control, default operator | t = 8.09 | 10.0 h | — | — |
| 20973341 | control, Yang operator | t = 13.7 | 20.8 h | — | stop ~00:40 Oct 1, t ≈ 15.6 |

There are no pressure-loop failures.

**Grain:**
- terminal settling at v = −1.969 = 0.922 w_s (Schiller–Naumann), constant since t ≈ 4;
- about 8% below the unbounded estimate, consistent with the hindrance from 8 d periodic
  neighbours;
- the path is exactly on axis (x = z = 4.0000);
- at t = 7.26 it is 5.3 d above the bottom wall, so expect contact around t ≈ 10.

**Grain-carried meltwater (first result).** Tracer integral below depth z, P1b minus the
default control, at t = 7.0:

| Deeper than (below ice) | P1b | Control | Excess (d³) | Grain volumes |
|---|---|---|---|---|
| 4 d | 4.195 | 3.601 | +0.593 | 1.13 |
| 6 d | 0.870 | 0.493 | +0.377 | 0.72 |
| 8 d | 0.216 | 0.041 | +0.175 | 0.33 |
| 12 d | 0.0398 | 6×10⁻⁵ | +0.0398 | 0.076 |

The grain carries about one grain volume of pure-meltwater equivalent out of the sublayer.
Its wake leaves a meltwater trail many diameters below the ~2 d zone affected by the band
kernel (§8.5).

**Health and budgets at t = 7.25:**
- salt ∫s ≤ 1.2×10⁻¹⁴;
- tracer − melt 1.4×10⁻⁷;
- heat drift 9.4×10⁻⁷;
- liquid salinity ≤ 1 and liquid tracer ≥ 0;
- net freezing continues slowly: −0.353 at t = 3, −0.422 d³ at t = 7.25 (default-kernel artifact);
- total wake ice 2.6×10⁻³ d³, about 0.5% of a grain volume.

### 8.7 Production segment complete — 2026-10-01

All three runs ended normally at the 23.5 h checkpoint-stop:

| Run | Reached | Cost |
|---|---|---|
| P1b (20977443) | t = 20.50 | 3,009 SU |
| Default control (20977444) | t = 18.32 | 3,010 SU |
| Yang control (20973341) | t = 15.43 | 3,009 SU |

The balance is now 62,279 SU. The full analysis is in
**[RESULTS_20977443.md](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf/RESULTS_20977443.md)**,
with figures F1–F5.

**Results:**
- **Settling:** terminal speed 14.35 mm/s (0.920 × Schiller–Naumann, hindered by the 8 d
  periodic neighbours), reaching 90% of it 2.65 d below the ice. Contact at t = 10.35
  (0.283 s), with no rebound.
- **Meltwater transport:** the grain carries about **one grain volume** of meltwater
  equivalent below 4 d (peak 1.26, 0.99 at t = 18.3).
  - The shallow part rebounds; the deep trail persists (0.47 below 8 d, 0.22 below 12 d).
  - The excess is centred 9.4 d below the ice.
- **Per-grain change to the mean stratification:** |ΔS| ≤ 0.003, |Δθ| ≤ 0.002.
- **Budgets:** roundoff-level salt (10⁻¹⁴), tracer 10⁻⁶, heat 10⁻⁹. Rigidity ≤ 0.65%
  throughout.

**Measured limitations:**
- the default band kernel gives net freezing (−0.34 d³) where the Yang reference melts
  (+0.34 d³);
- the near-ice bias reaches 3.4 d in θ by t = 15;
- the wake ice is ≤ 0.8% of a grain volume.

**Suggested next steps:**
1. repair the Yang operator at resolved sediment;
2. run a no-ice twin to separate ceiling hindrance from stratification drag;
3. run P2 at 32/d (about 9.4k SU at the measured 7.4k cell-steps per core-second);
4. run a concentration study.

### 8.8 Follow-up campaign — submitted 2026-10-01 (user: "proceed")

| Item | Job | Build | Purpose |
|---|---|---|---|
| Concentration variant `conc4` | **21006061** | v2 | Lateral period 4 d (one grain per 4 × 4 d), same physics; compared per area with the 8 d results |
| Thin-strip controls `strip24_def`, `strip32_def` | **21006078 / 21006079** | v2 | No-grain controls on a 4 × 4-cell strip |
| **P2**, 32 cells/d | **21006081** | v2 | Refined precursor paired with P1b (default operator) |
| dev3 (debug) | **21006068** | v2 / v3 | Repaired-Yang test and v2/v3 bit-identity check |
| No-ice twins `twin_lid_strat`, `twin_lid_uniform` | after dev3 | v3 | Separate ceiling hindrance from stratification drag |
| `p1b_yang` with control `strip24_yang` | after dev3 passes | v3 | Accurate near-ice and melt quantities |

Notes on each item:
- **Thin-strip controls:**
  - the no-grain control is horizontally uniform (its dissipation is 10⁻¹⁴), so a strip
    reproduces its horizontal mean at about 10⁻³ of the cost;
  - `strip24_def` is validated against the full control 20977444 before `strip32_def` is
    used as P2's control.
- **P2:**
  - timestep cap halved per B.3.5 (0.0025), run to t = 12 (settling, landing and trail);
  - run on 4 nodes, 512 ranks, `wholenode`: one node would need about 4 days of restart
    segments;
  - about 8.7k SU expected, 12.3k cap.
- **dev3:**
  - `Y_orig` (v2) is expected to diverge; `Y_fix` (v3) must keep liquid S within
    [s_i − 0.02, 1];
  - G2r (1-D Yang gate, no grain) must match round 2 to 10⁻⁶;
  - default-operator output must be bit-identical between v2 and v3;
  - criteria declared in `check_dev3.py`.
- **No-ice twins:** a no-slip lid at the same distance; one with the P1b sublayer, one with
  uniform T/S plus the same passive tracer.

**Code (build_v3, additive):**
- **Yang repair:** under `VOF_IBM` the Yang capacity is C_L + C_S, i.e. 1 − ice
  (`Conc_yang_capacity_field`, used in the implicit operator, the explicit diffusion, the
  inner product and the dilution source).
  - Rationale: sediment is not ice, so the support cells no longer get a 1/δ ≈ 10⁶
    coefficient.
  - It is identical to C_L when there is no grain.
- **New key `[conc] sublayer_y0`:** an optional anchor for the init 36–38 sublayer profiles,
  default inert.

**Update 2026-10-01, 16:55.**
- **First strip attempt failed:** 21006078/79 stalled in HYPRE (residual ~0.75 on a 4-cell
  strip with buoyancy on).
- **Resubmission:** 21006850/51, with buoyancy **off**. The full no-grain control is
  quiescent, so its buoyancy is pure hydrostatics.
- **`strip24_def` validated** against the full control 20977444 at t = 1, 3, 7, 12 and 18:
  per-area tracer, ice, salt and T integrals, plus mean T/S, agree to **≤ 1.6×10⁻¹⁰**
  (mean T ≤ 8×10⁻¹⁰), relative 10⁻¹⁰–10⁻⁸.
- **`strip32_def`** (P2's control) finished to t = 12 in 2 minutes.
- Lesson: a no-grain control in this configuration costs minutes, not a 3k-SU full box.
- Still queued: dev3 (debug) and P2. `conc4` is running (t = 0.49 at 17 min, one pressure
  pass).

### 8.9 Follow-up results so far — 2026-10-02

**dev3 (21006068) passes all three pre-declared criteria:**
- **R1:** the repaired Yang operator (build_v3) keeps liquid S within [0.6276, 1.0000] to
  t = 1.2 with no pressure failure. On v2 the same deck diverged (S ±11 at t = 0.6) and
  failed at t = 0.734, reproducing P1b-Yang.
- **R2:** the 1-D Yang gate without a grain gives λ = 0.019412720023085504, identical to
  round 2.
- **R3:** default-operator output is bit-identical between v2 and v3.

So the Yang operator is usable with moving grains in build_v3.

**conc4 (21006061, 982 SU):** per-grain meltwater transport falls by about 40% when the grain
spacing goes from 8 d to 4 d, and terminal speed by 18%. Transport is sub-linear in
concentration ([RESULTS](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf/RESULTS_20977443.md)
§2b, F6).

**Submitted 2026-10-02 on build_v3:**
- no-ice twins **21019097** (lid with stratification) and **21019098** (lid with uniform
  water);
- **P1b-Yang (repaired) 21019100** with strip control **21019101**.

**Still queued:** P2 (21006081, 4 nodes).

## 9. Reproduce

```bash
cd PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf
/apps/anvil/external/apps/anaconda/2025.06/bin/python3 prepare.py        # all decks + case.json
python analyze_gate.py <gate_smoke_dir>/gate                              # gate verdict
python check_smoke.py <gate_smoke_dir>/smoke --case case.json             # smoke checks
# stage a copy under ECCO_StageB3/P1_Attempts/P1b_shelf/, then:
sbatch job.sh ; sbatch --export=ALL,CASE=ctrl job.sh
```

## References

- A. Jenkins, K. W. Nicholls, H. F. J. Corr (2010), [Observation and parameterization of ablation at the base of Ronne Ice Shelf](https://journals.ametsoc.org/view/journals/phoc/40/10/2010jpo4317.1.xml), *JPO* 40.
- T. Keitzl, J. P. Mellado, D. Notz (2016), [Reconciling estimates of the ratio of heat and salt fluxes at the ice–ocean interface](https://arxiv.org/abs/1606.03004).
- D. Shean et al. (2019), [Pine Island ice-shelf basal melt rates](https://tc.copernicus.org/articles/13/2633/2019/tc-13-2633-2019.html).
- E. Pierce, I. Overeem, B. Hasholt (2025), [Sediment transport by Greenland's icebergs](https://pmc.ncbi.nlm.nih.gov/articles/PMC12859035/).
- Srdić-Mitrović, Mohamed & Fernando (1999) and Abaid et al. (2004): particles crossing density interfaces.
- [TEOS-10 GSW density](https://teos-10.github.io/GSW-Python/density.html).
