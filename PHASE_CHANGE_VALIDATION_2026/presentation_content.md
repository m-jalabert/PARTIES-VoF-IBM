# Presentation content — PARTIES phase-change implementation
Source: main_PHASE.tex (validation report) + Boundary.h flags + NASA-JPL/YANG_ECCO_IMPLEMENTATION_ROADMAP.md
Each block below = one slide. Copy into your template.

---

## Slide 1 — Title

**A diffuse-interface phase-change model for melting in salt-stratified convection**
Implementation and validation in PARTIES

- M. Jalabert — UCSB, Mechanical Engineering
- Three benchmarks: Stefan (verification) → melting Rayleigh–Bénard (validation) → Yang et al. 2023 (application)

---

## Slide 2 — Motivation

- Goal: simulate sediment release from ice melting into salt-stratified water (motivates an ECCO ocean-science application)
- This requires **three physics coupled at once**: a moving phase boundary, double-diffusive convection, and a nonlinear (cold-water) equation of state
- None of this existed in the base PARTIES multiphase solver — it is a new capability, added as **four independent, flag-gated features**
- Non-invasiveness constraint: with all four flags off, the code is bit-for-bit identical to the existing multiphase solver

---

## Slide 3 — Model overview: one-fluid diffuse-interface formulation

- Boussinesq, one-fluid: ice and water share ρ₀ → velocity stays solenoidal, no density jump in the pressure projection
- Phase indicator F = liquid fraction (Cahn–Hilliard field), solid fraction φₛ = 1 − F
- Nondimensionalization: length by H, velocity by free-fall U = √(gΔρH/ρ₀), θ = (T−Tᵢ)/ΔT, s = S/Sₘ
- Four new physics blocks, each behind its own compile flag in `Boundary.h`:

| Flag | Physics added |
|---|---|
| `PHASE_CHANGE` | Stefan melting source in the Cahn–Hilliard equation + matching latent-heat sink in temperature |
| `CONC_VOF_PHASEWEIGHTED` | F-weighted scalar diffusivity — salt confined to the liquid, heat diffuses through both phases |
| `EOS_NONLINEAR` | Roquet quadratic buoyancy b(θ,s) — density maximum inside the temperature range |
| `ICE_PENALIZATION` | Darcy/Brinkman damping that holds the unmelted solid rigid |

- Each benchmark below turns on a **different subset**, building up to the full stack

---

## Slide 4 — Governing equations (as solved)

Continuity + momentum:

  ∇·u = 0
  ∂u/∂t + ∇·(uu) = −∇p + √(Pr/Ra) ∇·[μᵣ(F)(∇u+∇uᵀ)] − b(θ,s) ŷ − (φₛ/τ) u

Interface (conservative Cahn–Hilliard + melt source):

  ∂F/∂t + ∇·(uF) = (1/Pe_CH) ∇²ψ + m,     ψ = W′(F) − Cn²∇²F
  m = (St/(Pe_T ε)) [θ − θ_L(s)] w(F),     w(F) = F(1−F)/(√2 Cn),     θ_L(s) = T_m − Λ*s

Temperature / salinity:

  ∂θ/∂t + ∇·(uθ) = (1/Pe_T) ∇·[κᵣ(F)∇θ] − m/St
  ∂s/∂t + ∇·(us) = (1/Pe_S) ∇·(F∇s)

Discrete enthalpy (the single most useful diagnostic, tracked in every case):

  E = ∫(θ + F/St) dV   — conserved by construction because m appears identically in the F and θ equations

---

## Slide 5 — Benchmark roadmap

| Benchmark | Flags exercised | Role |
|---|---|---|
| 1. Binary Stefan (1D) | `PHASE_CHANGE`, `CONC_VOF_PHASEWEIGHTED` | Verification against an exact similarity solution |
| 2. Melting Rayleigh–Bénard | `PHASE_CHANGE`, `ICE_PENALIZATION` | Validation against published DNS (Favier et al.) |
| 3. Yang et al. salt-stratified melting | `PHASE_CHANGE`, `CONC_VOF_PHASEWEIGHTED`, `ICE_PENALIZATION`, `EOS_NONLINEAR` (full stack) | Application-relevant validation, double-diffusive convection |

Each benchmark adds one more piece of physics on top of the last — by benchmark 3 every flag is active simultaneously.

---

## Slide 6 — Benchmark 1: Binary Stefan problem — setup

**What's tested:** the Stefan source and latent sink in isolation (`PHASE_CHANGE`), coupled to the liquidus depression and phase-weighted salt confinement (`CONC_VOF_PHASEWEIGHTED`). No convection, no buoyancy — the sharpest available test of the melting physics itself.

- Planar interface separates saline liquid from fresh ice; melting temperature depressed locally by salinity
- Exact self-similar solution: front position X(t) = X₀ + 2λ√(t/Pe_T), λ solves a transcendental equation coupling the thermal and solutal similarity variables (λ, λₛ = λ√Le)
- For the chosen inputs: λ = 0.21713827, sᵢ = 1.16402×10⁻³, θᵢ = −5.82010×10⁻⁴
- Flow is exactly quiescent (max|u| = 0) — isolates the source terms from any advective error

---

## Slide 7 — Benchmark 1: Binary Stefan — inputs

| Parameter | Value |
|---|---|
| Domain Lx×Ly | 1.5 × 7.8125×10⁻³ |
| Grid Nx×Ny×Nz | 768×4×1 |
| Thermal Péclet Pe_T | 100 |
| Solutal Péclet Pe_S | 10⁴ |
| Lewis number Le | 100 |
| Stefan number St | 0.5 |
| Liquidus slope Λ* | 0.5 |
| Cahn number Cn = band scale ε | 1.46484×10⁻³ |
| CH Péclet Pe_CH | 614.4 |
| Far-field θ∞, s∞, ice θ_ice | 1, 1, 0 |
| End time | 0.5 |

---

## Slide 8 — Benchmark 1: Binary Stefan — results

| Quantity | Computed | Exact | Rel. error |
|---|---|---|---|
| Similarity constant λ | 0.2184496 | 0.2171383 | 6.0×10⁻³ |
| Front-fit RMS residual | 6.3×10⁻⁶ | — | — |
| Final salinity L∞ error | 1.44×10⁻² | — | — |
| Salt conservation drift | 6.0×10⁻¹⁵ | 0 | — |
| Enthalpy drift | 5.0×10⁻¹⁰ | 0 | — |
| Max velocity | 0 | 0 | — |

- Front follows √t from the first output; compositional boundary layer (~7×10⁻³ thick, 10× thinner than the thermal layer) captured without oscillation
- Salinity drops to sᵢ ≈ 1.2×10⁻³ at the front — released meltwater is fresh
- **Takeaway:** the Stefan source, latent sink, and liquidus coupling reproduce the exact solution to <1%, with machine-precision conservation

---

## Slide 9 — Benchmark 2: Melting Rayleigh–Bénard — setup

**What's tested:** `ICE_PENALIZATION` (Brinkman rigidity) in a genuinely convecting flow, on top of `PHASE_CHANGE`. First case where the phase boundary interacts with buoyant convection instead of sitting in a quiescent layer.

- Favier, Purseed & Duchemin benchmark: solid layer melts from below into a convecting fluid
- Two regimes: an early **conduction-dominated** phase with a known constant (½+St_alt)Pe·ḣ, and a later **convective plateau** where the effective Nusselt number follows Ra_e^(1/3)
- Freshwater only (Λ* = 0) — isolates thermally driven melting from the salt physics tested in Benchmark 1
- Runs match the reference's case D (Ra=10⁷, aspect ratio 6, θ_M=0.05), at 3× finer resolution in x

---

## Slide 10 — Benchmark 2: Melting Rayleigh–Bénard — inputs

| Parameter | Value |
|---|---|
| Domain Lx×Ly | 6 × 1 |
| Grid | 3072×512×1 |
| Reynolds Re = √(Ra/Pr) | 3162.28 |
| Thermal Péclet Pe_T | 3162.28 |
| Stefan number St | 10 (St_alt = 0.1) |
| Melting temperature T_m | 0.05 |
| Liquidus slope Λ* | 0 (no salt) |
| Cahn number Cn = ε | 1.46484×10⁻³ |
| CH Péclet Pe_CH | 614.4 |
| Darcy time τ | 10⁻⁴ |
| Initial solid thickness | 0.05 |
| End time | 120 |

---

## Slide 11 — Benchmark 2: Melting Rayleigh–Bénard — results

- **Diffusive phase:** measured (½+St_alt)Pe·ḣ = 23.34 vs analytical 23.14 → **0.9% agreement**; max deviation from the exact 1D Stefan solution over that phase: 1.7%
- **Convective plateau:** γ_eff = Nu/Ra_e^(1/3) = 0.0983 vs published 0.115 (ratio 0.854); interface energy balance closes to 0.996
- **Stefan-number scaling:** measured ḣ(St_alt=0.1)/ḣ(St_alt=1) = 3.01 vs expected 2.5
- **Direct DNS comparison** (developed melting velocity vs St, Favier fig. 7b): −4.8% (St_alt=0.1), −15.3% (St_alt=1) on a quantity spanning two orders of magnitude
- **Heat transport law:** both melting runs and an independent flat-wall Ra=10⁶ point collapse onto a common Ra_e^(1/3) branch — prefactor 0.093–0.097 vs reference 0.115 (16–19% low), but the **scaling exponent is reproduced** and melting doesn't perturb the transport law
- **Morphology:** front flat while conducting, develops cusps once convection sets in (t≈20), each cusp locked to a downwelling; 14 extrema across Lx=6 at t=70, mean spacing 0.86 ≈ fluid depth
- Flat-wall control: Nu = 8.27±0.26 at Ra=10⁶, confirming no convective-transport deficit underlies the melting cases

---

## Slide 12 — Benchmark 3: Yang et al. (2023) salt-stratified lateral melting — setup

**What's tested:** the full flag stack together — `EOS_NONLINEAR` added on top of `PHASE_CHANGE` + `CONC_VOF_PHASEWEIGHTED` + `ICE_PENALIZATION`. This is the application-relevant case: double-diffusive convection with a nonlinear, cold-water equation of state driving lateral ice melt.

- Vertical ice block (thickness 0.1H) melts laterally into stably salt-stratified water, closed insulated square cavity
- Two cases: freshwater (Sₘ=0, isolates thermal melting) and salt-stratified reference (Sₘ=5 g/kg, ΔSᵥ=5 g/kg)
- Buoyancy: quadratic Roquet EOS b(θ,s) = −βT|θ−T_md,0−σs|² + βS s — places a density maximum inside the temperature range, essential for cold-water dynamics
- Key comparison-protocol issue: the **reference uses a dual-grid scheme** (velocity/temperature on 288², salinity/phase field refined 5× to 1440²); PARTIES is single-grid, so no single run matches both simultaneously

---

## Slide 13 — Benchmark 3: Yang — inputs

| Parameter | Value |
|---|---|
| Domain, aspect ratio | H×H, Γ=1 |
| Initial ice thickness | 0.1H at x≥0.9H |
| Thermal Rayleigh Ra_T | 10⁷ |
| Prandtl Pr | 10 |
| Schmidt Sc (Le=100) | 10³ |
| Reynolds Re = √(Ra_T/Pr) | 10³ |
| Pe_T, Pe_S | 10⁴, 10⁶ |
| Stefan number St | 0.25 (St_alt=4) |
| Liquidus slope Λ* † | 0.014 |
| EOS: βT, βS †, T_md,0, σ † | 1, 1.75, 0.2, −0.0625 |
| Cahn number Cn | 0.75Δx |
| CH Péclet Pe_CH | 0.9/Cn |
| Darcy time τ | 10⁻⁴ |
| BCs | no-slip, no-flux, all walls |
| End time | 200 |

† S_m-dependent, vanish in the freshwater case

---

## Slide 14 — Benchmark 3: Yang — resolution protocol finding

- Measured resolution dependence directly, N=216→1440, identical physics/binary — only mesh + 4 mesh-tied numbers vary

| N | 216 | 288 | 360 | 512 | 720 | 1440 |
|---|---|---|---|---|---|---|
| Freshwater t½ | 116.31 | 102.36 | 94.47 | 85.59 | 78.89 | 60.29 |
| Salt-stratified t½ | — | 188.39 | — | 178.28 | — | 179.85 |

- **Freshwater t½ varies 70%** across the ladder — convection-driven, resolution-limited heat delivery
- **Salt-stratified t½ varies only 5.4%**, non-monotonically — the vertical salinity gradient suppresses convection enough that 288² already resolves it
- Consequence: comparing at the reference's *finest* grid (1440²) vs its *coarsest* (288², which is what actually carries velocity/temperature) gives very different apparent agreement — the protocol dominates the comparison unless resolution dependence is measured first

---

## Slide 15 — Benchmark 3: Yang — results at matched resolution

Comparison at the reference's effective resolution (288², since it carries u,θ there in both cases):

| Quantity | PARTIES 288² | Reference | Difference |
|---|---|---|---|
| Freshwater t½ | 102.36 | 107.42 | −4.7% |
| Freshwater V(200)/V₀ | 0.2447 | 0.2483 | −1.4% |
| Salt-stratified t½ | 188.39 | 216±11 | **−12.8%** |
| Melt-rate ratio f̄₅/f̄₀ | 0.543 | 0.497 | +9.3% |

For contrast, at protocol-mismatched 1440²: freshwater t½ = −43.9%, melt-rate ratio = −32.6%

- Field structure: 3 stratification-driven horizontal layers, each with its own convective cell against the ice — matches reference figs. 1–2
- Interfacial velocity: freshwater downward jet peaks at |v|=0.080 vs reported ≈0.085; stratification weakens it 4×, to 0.021
- Front shape: freshwater smooth (0.0056H peak-to-peak), stratified scalloped at layer spacing (0.019H, 3.4× larger)
- **Open item:** the −12.8% salt-stratified residual is *not* a resolution artifact (that case is converged) — it is real and reported as unexplained (2.6σ against a reference value that is itself derived, ±4.9%)

---

## Slide 16 — Next step: trapped/resolved sediment particle (Stage B)

**Motivation:** the original driver for this work — sediment release from melting ice — requires a resolved particle held rigidly inside the ice and released when the surrounding ice melts. This is the natural extension once the Stage-A melting physics (Slides 6–15) is validated.

- Combines `VOF_IBM` + `LAG_PARTICLE_RESOLVED` with the Stage-A phase-change flags — gated so Stage-A behavior is provably unchanged (bit-identical regression passed)
- **Hold-in-ice:** the grain lives in the *mobile* particle list from t=0 with its motion locked (U=Ω=0), not moved between fixed/mobile lists — avoids corrupting the HDF5 output layout and lets it accumulate a converged hydrodynamic force history before release
- **Release criterion:** shell-averaged liquid fraction φ_liq around the grain (not a centroid probe — the grain's interior always reads "solid"); released once φ_liq ≥ F_release, then ramped to free motion
- **Meltwater tracer:** a 3rd, non-buoyant scalar injected at the interface at the local melt rate — makes entrainment tracking a built-in conservation check (∫C_mw should track ∫F exactly)
- **Found during shakeout:** the VOF admissibility clip that keeps F∈[0,1] silently destroys liquid mass where it overlaps the sediment's solid indicator (sediment cells are exempt from the usual conservative redistribution) — a genuine budget leak, not a diagnostic artifact, growing with grain surface area
- **Status:** 2D release shakeout passed after fixing three bugs (predictor–corrector lock leak, rank-local diagnostic, restart persistence); the IBM/VOF mass leak is addressed by extending local conservation into the diffuse band (same mechanism as the existing band mass-restore fix)
- **Remaining:** 3D debugging pass (single particle, reduced domain) → full ECCO production campaign (Galileo-number sweep, Sc=7/70 extrapolated to oceanic Le≈100) → offline mixing diagnostics (κ_eff, background-PE mixing efficiency), already implemented and smoke-tested

---

## Slide 17 — Summary

- **Stefan (1D):** exact-solution verification of the Stefan source + liquidus coupling — front constant to 0.6%, salt/enthalpy conserved to machine precision
- **Melting Rayleigh–Bénard:** convective coupling validated against published DNS — diffusive phase to 0.9%, heat-transport scaling exponent reproduced, no convective-transport deficit
- **Yang et al.:** full coupled stack (nonlinear EOS + salt confinement + convection + ice rigidity) — freshwater matched to 4.7%/1.4%, salt-stratified melt rate to 9.3%, with a genuine 12.8% residual reported openly rather than hidden by a favorable but mismatched resolution comparison
- **Methodological point:** apparent agreement with a reference depends heavily on matching its resolution protocol — always measure your own resolution dependence before attributing a gap to the physical model
- **Next:** extend to a resolved sediment particle trapped in and released from melting ice, the ECCO-motivated application
