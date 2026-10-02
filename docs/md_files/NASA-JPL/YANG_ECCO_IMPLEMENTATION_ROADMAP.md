# Implementation Roadmap — Yang Melting Benchmark → ECCO Sediment-from-Ice DNS

> **B.2 CLOSED, 2026-09-08.** All debugging checks pass, re-verified under
> a solver fix in job **20492853**. The [closure report](YANG_ECCO_B2_CLOSURE_REPORT.md)
> supersedes the historical handover and old Darcy/confinement recommendations.
> Actual-ice transport, conservative moving-mask remapping, and periodic particle
> geometry are fixed and promoted behind default-off switches. **Use tau=1e-3
> with the hard threshold in the validated configuration**, not the old 0.391.
> Corrected settling refinement, budgets, MPI remap, restart, profiles, clean-build
> reproducibility, and final Stage-A regression pass. The old terminal-velocity/
> Faxén claim is withdrawn; the short column establishes transient refinement.
>
> **The matched-ramp timestep comparison 20487690 was run and FAILS.** It found a
> genuine solver defect — the release ramp re-damped the accumulated velocity each
> RK stage, so the post-ramp velocity scaled as `sqrt(dt)` with no convergent limit
> (fixed; `rampfix.patch`) — and, more importantly, a limitation refinement cannot
> remove: **the grain's escape from the ice is an ill-conditioned near-cancellation
> and its timing does not converge.** The IBM reaction supports 97.9-99.9% of the
> buoyant weight during the creep, so the escape time drifts (5.75 -> 6.63 -> 7.42
> as dt halves, observed order 0.15) and coarse timesteps release the grain
> systematically early. This is not chaos: a 64-vs-32-rank twin is bit-identical,
> and pre-release fields converge at order ~0.9.
>
> **Constraint on B.3:** do not report or depend on absolute release time, escape
> duration, or release-to-contact timing; the bias is one-signed, so an ensemble
> does not remove it. Ice-free settling velocity IS converged to 0.28% at
> max_dt=0.01. Separately, the hard ice-penalization threshold costs ~8x in
> settling dt-convergence (3.29% vs 0.41% with it disabled), but disabling it
> delays release by 10.3 time units and worsens ice rigidity — a documented
> trade-off, not a recommended setting.
>
> **B.3 RESULT, audited 2026-09-27:** P1 feasibility job **20906327 COMPLETED**
> for **182.76 SU**, reaching t=10.003605 (0.61994 physical seconds). Conservation
> passes, but there is **no release**, slight **net freezing**, and the relative
> ice-rigidity criterion **fails (3.11% versus 1%)**. This completes the bounded
> feasibility segment, **not full P1 acceptance**. See §B.3 and the
> [result report](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/RESULTS_20906327.md).
> No continuation, P2 or production was submitted during this audit.
>
> **B.3 REDESIGN, 2026-09-28 (controlling): P1a is superseded by P1b.** P1a's net
> freezing was a diffuse-band artifact rather than the physics of its inputs. The
> exact binary-Stefan solution melts (λ = 0.0112). The band read the volume-averaged
> salinity at a liquidus slope of about 1, so it pinned at θ ≈ 0.44 instead of 0.05.
> Waiting for a trapped grain to melt out would cost about 10⁵–10⁶ SU for any
> realistic polar forcing. **Release is now the initial condition.** P1b places a
> free 0.2 mm quartz grain beneath an Antarctic ice-shelf base (700 dbar, mixed layer
> 0 °C / 34.85 g/kg, u* = 0.5 cm/s), inside the three-equation diffusive meltwater
> sublayer (interface 21.30 g/kg, −1.670 °C, 35 m/yr). Melting is ON, with a matched
> ice-only control.
>
> **Update 2026-09-30.** Two pre-flight rounds replaced the band treatment with the
> **Yang finite-interface salt operator**, and P1b (**job 20973340**) plus its control
> (**job 20973341**) were submitted.
> - The 1-D gate rejected both default-operator variants (λ +94% and −9%). The Yang
>   operator is within 0.2%, but its salt drift is 1.8×10⁻⁵ against a pre-declared 10⁻⁸;
>   that deviation is recorded.
> - `liquid_referenced_salinity = 1` also stalls the pressure solve with a resolved grain,
>   because it divides by C_L ≈ 0 inside the sediment.
>
> **Update 2026-09-30, 11:30.** P1b-Yang (20973340) **failed at t = 0.738**: the Yang salt
> field diverges at a moving resolved grain, because its 1/(F+δ) factor reaches ~5×10⁵ at the
> C_L = 0 support cells and Crank–Nicolson does not damp it.
> - P1b was resubmitted with the **default operator**, which is B.2-validated for moving
>   grains; its known cost is band-local (λ −9%, wrong interface salinity). New jobs:
>   **20977443** plus control **20977444**.
> - The Yang control 20973341 continues, to measure how much the band-kernel choice changes
>   the sublayer.
>
> **P1b RESULT, 2026-10-01.** The production segment completed: P1b reached t = 20.5 and the
> controls t = 18.3 (default) and 15.4 (Yang), for about 9.0k SU in total; 62,279 SU remain.
> - **Settling:** terminal 14.35 mm/s (0.92 × Schiller–Naumann, box-hindered), landing at
>   0.28 s.
> - **Meltwater transport:** the grain drags **about one grain volume of meltwater** below
>   4 d, and a persistent trail of 0.22 grain volumes below 12 d.
> - **Budgets:** roundoff-level.
> - **Limitations:** near-ice results (≲ 3.4 d), the melt rate and its sign are
>   band-kernel-biased.
>
> See [RESULTS_20977443.md](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1b_shelf/RESULTS_20977443.md).
> P2 and the Yang-at-sediment repair are proposed, not submitted.
>
> The dedicated record is
> [YANG_ECCO_B3_P1b_SHELF.md](YANG_ECCO_B3_P1b_SHELF.md). `testcases/ECCO_TESTS/`
> was reorganized by stage the same day (see its `README.md`).

**ABSOLUTE GOAL (controlling; restated 2026-07-27).** Reproduce the Yang et al. (2023) benchmark
**with PARTIES**. Running the authors' own solver is **not** a deliverable and is **not** an
acceptable substitute: it would validate AFiD-MuRPhFi, not PARTIES, and it does nothing for
Stage B, which is a PARTIES-only sediment-from-ice study. Every simulation reported for this
benchmark is run with PARTIES.

> **Withdrawn recommendation.** Between 2026-07-18 and 2026-07-21 this roadmap (§A.4 closure
> block), [YANG_SALT_OPERATOR_GATE_REPORT.md](YANG_SALT_OPERATOR_GATE_REPORT.md) §7, and
> [YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md) §7.5 recommended "reproduce Yang
> with AFiD-MuRPhFi" as the allocation-minimum route, and demoted a PARTIES result to a
> conditional second option requiring a full Allen–Cahn/Hester port. **That recommendation is
> withdrawn as off-goal.** The evidence in those reports stands; only the recommended action
> changes. AFiD remains usable as a *source-code reference* for algorithmic details and its
> published 1-D validation remains usable as a *target* for a PARTIES gate — but no AFiD run is
> part of this project. The replacement closure plan is **§A.5**.
>
> The withdrawal is also justified on the merits, not only on scope: nothing in the accumulated
> evidence shows that PARTIES' conservative Cahn–Hilliard formulation is the defect (see §A.5.0),
> and the one head-to-head test that exists — the Le=100 salt-operator gate against an exact
> binary-Stefan similarity solution — was **won by the PARTIES formulation**, not by the mapped
> Yang operator. Porting Allen–Cahn is therefore a hypothesis to be tested last, not a
> prerequisite.

**Scope of this document.** A lean, two-stage plan:

1. **Stage A — Yang benchmark.** All code modifications (written 3D), a *debugging* pass in a
   tiny reduced domain (code shake-out only), and **four preliminary 1D validations** — the
   1D Stefan (Neumann) gate (§A.3, **passed 2026-07-02**), the 1D binary/salty Stefan gate
   (§A.3b, **passed 2026-07-03**), the nonlinear-EOS reconstruction gate (§A.3c,
   **passed 2026-07-05**), and the ICE-penalization sweep (§A.3d,
   **passed 2026-07-05**) — then the **single full 2-D Cartesian
   production run** that reproduces Yang et al., *J. Fluid Mech.* **969**
   (2023), R2 — "Ice melting in salty water: layering and non-monotonic dependence on the mean
   salinity"
   ([DOI 10.1017/jfm.2023.582](https://doi.org/10.1017/jfm.2023.582)·
   [arXiv:2302.02357](https://arxiv.org/abs/2302.02357)).
2. **Stage B — ECCO study.** The additional code updates (resolved sediment particle held in
   ice, interface-triggered release, meltwater tracer, mixing diagnostics), a debugging pass in a
   **reduced affordable 3D domain**, then the **full 3D ECCO production inputs**.

**Conventions.** Code storage is 3D, but the Stage-A production benchmark uses `TWOD_CARTESIAN`.
PARTIES gravity is `−y` ⇒ **vertical = y**. Coordinate map
from Yang (their vertical `z`, horizontal `x`, depth `y`): **Yang z → PARTIES y**, **Yang x →
PARTIES x** (ice slab at `x ≥ 0.9H`), **Yang depth y → PARTIES z** (storage-only slab in the
2-D benchmark).
Nondimensionalization is thermal free-fall: length `H`, `U = √(g·(Cb ΔT²/2ρ₀)·H)`,
`θ = (T−Ti)/ΔT ∈ [0,1]`, `s = S/Sm`.

**Deliberate omissions (per request).** No Rayleigh–Bénard, no lock exchange, no
penetrative-convection test, no incremental V1–V11 ladder. Stage A keeps **exactly four**
preliminary physics validations, all cheap 1D gates with exact references:
two quiescent Stefan gates — the **1D Stefan (Neumann) test** (§A.3, passed 2026-07-02), which
certifies the phase-change kernel alone (G.2 source + G.3 condition + G.4 latent sink), and the
**1D binary (salty) Stefan test** (§A.3b, passed 2026-07-03), which additionally certifies the
salt–liquidus coupling (G.3 liquidus, G.5 salt masking, solid-side conduction through
`κ_r(F)`) — and two **steady unidirectional-flow gates** (added 2026-07-05) for the two
flow-coupled modules the quiescent gates cannot reach: the **EOS reconstruction test** (§A.3c),
which measures the Roquet buoyancy `b(θ,s)` pointwise from the curvature of a steady velocity
profile (G.6), and the **Brinkman ice-penalization test** (§A.3d), which measures the Darcy
damping magnitude, its `δ = √(ν·τ)` boundary-layer, and the `√τ` slip scaling (G.1). Together
the four gates cover every new Stage-A module before the expensive Yang run; everything else in
Stage A is the full 2-D Yang case. The only physics run in Stage B is the ECCO production
case. "Debugging" below means short, tiny-grid runs whose sole purpose is to shake out
compile/run/MPI/conservation bugs — not to validate physics; the 1D gates (§A.3–§A.3d), by
contrast, **are** quantitative physics validations.

**Parallelization requirement (applies to every code modification below).** PARTIES runs
distributed-memory (MPI pencil decomposition + HYPRE). **All new code must be multiprocessor-safe**
— it has to give bit-for-bit the same result on any process count and add no serial bottleneck.
Concretely, every new term must:
- operate only on the local sub-domain and obtain off-rank neighbours through the existing **ghost
  /halo exchange** (≥ 1 ghost layer; the `∇F` stencil of the melt source in G.2–G.3 needs a fresh
  halo update **before** it is evaluated);
- use **collective MPI reductions** for any global quantity (ice volume `V(t)`, heat/salt/tracer
  budgets, melt-rate normalization, `κ_eff`, and especially the **BPE density sort** in B.1.6, which
  must be a parallel sort/gather — never a serial gather to rank 0);
- respect the pencil decomposition in the new initial conditions (each rank fills only its slab; no
  global allocation), and keep new `[eos]`/`[phase_change]` inputs broadcast to all ranks;
- preserve HYPRE/CG compatibility — the Darcy penalization (G.1) and phase-weighted diffusivity
  (G.4–G.5) enter the **existing** implicit operators as spatially-varying coefficients, so the
  matrix assembly stays domain-decomposed.

Each modification below is annotated with **[MPI]** where it touches halos or global reductions.

**Non-invasiveness requirement (applies to every code modification below).** The new physics must
**not change the behaviour of the existing code**. Every addition is gated behind a **new
compile-time flag in `PARTIES/src/Include/Boundary.h`** (`EOS_NONLINEAR`, `CONC_VOF_PHASEWEIGHTED`,
`PHASE_CHANGE`, `ICE_PENALIZATION`, and the Stage-B particle flags), all **undefined by default**,
so that:
- with the new flags off, the code compiles and runs **bit-for-bit identically** to the current
  branch — every existing testcase still passes;
- new terms are wrapped in `#ifdef … #endif` blocks (or behind a runtime branch that is unreachable
  when the flag is off) and add **zero cost** to runs that don't enable them;
- new input keys (`[eos]`, `[phase_change]`, `darcy_tau`, `kappa_ice_ratio`, `F_release`) live in
  the `src/IO/default.inp` X-macro registry (so they are always *parsed*, with inert defaults) but
  are **used** only under their flag, and the new VOF/`Conc` `init_type` cases are **additive** —
  existing `init_type` values are untouched;
- no existing function signature, struct layout, or default input value is modified; new fields/
  members are appended, not inserted, and guarded by their flag.

This keeps the melting/ECCO development fully separable from the production multiphase solver until
each piece is validated.

---

## Governing equations (one-fluid / phase-field form, free-fall nondimensional)

All equations are written exactly as PARTIES solves them, in the thermal free-fall
nondimensionalization of the Conventions block (`θ = (T−Ti)/ΔT`, `s = S/Sm`, length `H`, velocity
`U`). Phase indicator `F` = liquid (water) volume fraction; the ice fraction is `φ_s = 1 − F`.
Each block ends with **→ implemented by** linking the new term to its code modification and inputs.
This is the Stage-A (Yang) system; Stage B adds the sediment-particle terms in §B.

### G.1 Continuity and momentum (Boussinesq one-fluid, with ice penalization)

The ice and water share the reference density ρ₀ (ρ₂/ρ₁ = 1), so the field stays solenoidal and the
existing pressure solver is untouched:

```
∇ · u = 0

∂u/∂t + ∇·(u u) = −∇p + √(Pr/RaT) ∇·[ μ_r(F) (∇u + ∇uᵀ) ]  −  b(θ,s) ŷ  −  (φ_s/τ_p) u
```

- `√(Pr/RaT) = 1/Re = 1/1000` (reference): viscous coefficient in free-fall units.
- `μ_r(F) = μ(F)/μ_w`: mixture viscosity, with `μ_ice/μ_w ~ 10²` as a rigidity backup.
- `− b(θ,s) ŷ`: buoyancy body force (ŷ = up); `b` is the nondim density anomaly `ρ'/Δρ_T` of G.6,
  so denser fluid (`b > 0`) accelerates downward.
- `− (φ_s/τ_p) u`: **Darcy / volume penalization** active only in ice (`φ_s = 1−F`), keeping the
  unmelted ice rigid.
- **Numerical treatment:** do **not** use a purely explicit RK source for production unless
  `τ_p` is intentionally weak (`τ_p` comparable to or larger than the RK substep). The production
  implementation should treat the drag as a local implicit Brinkman term in the face-centered
  velocity Helmholtz solve:

  ```
  [ ρ_f/(α_k Δt) + 2 ρ_f φ_{s,f}/τ_p ] u_f
    − (1/Re) ∇·( μ ∇u_f ) = RHS_nonpenalized
  ```

  The added coefficient is positive and diagonal, so the CG velocity operator remains SPD and the
  pressure Poisson equation is unchanged. Add the same diagonal term to the Jacobi preconditioner.
  When `LAG_PARTICLE_RESOLVED` is enabled in Stage B, use the same local diagonal in the explicit
  IBM predictor velocity so the IBM force is estimated against the damped ice-fluid state.
- **→ implemented by** A.1.5 (`ICE_PENALIZATION`, input `darcy_tau`) for the Darcy term and the
  `μ2/μ1` ratio; buoyancy by A.1.2.

### G.2 Interface evolution (Cahn–Hilliard) with Stefan melting source

The existing conservative CH transport gains a non-conservative melting source:

```
∂F/∂t + ∇·(u F) = (1/Pe_CH) ∇·( M ∇ψ )  +  V_Γ |∇F|
```

- `(1/Pe_CH)∇·(M∇ψ)`: existing CH right-hand side (mobility `M`, chemical potential `ψ`, interface
  thickness set by `Cn`). No surface-tension force is applied to the momentum equation — `ψ` only
  regularizes the tanh interface.
- `V_Γ |∇F|`: melting source localized on the interface; `V_Γ > 0` ⇒ the front recedes into the ice.
- **→ implemented by** A.1.4 (`PHASE_CHANGE`): new `VOF_DIFFUSE_compute_melt_rate` adds `V_Γ|∇F|`
  to `rhs_explicit` in `VOF_DIFFUSE_step`.

### G.3 Stefan condition and liquidus (superheat form, as implemented)

The energy jump at the interface supplies `V_Γ`. Dimensionally,
`ρ₀ L V_Γ = k_l ∂T/∂n|_liq − k_s ∂T/∂n|_sol` (with `St = cp ΔT/L`, the PARTIES convention;
`St_Yang = 1/St`). **The implemented form is the superheat-driven phase-field coupling of
Hester et al. (2020) / Yang et al. (2023):**

```
V_Γ = (St / Pe_T) ( θ − θ_L(s) ) / ε ,     ε = melt_band_eps  (≈ Cn)

θ_L(s) = T_melt − Λ* s ,     Λ* = liquidus_slope = m_freeze · Sm / ΔT
```

evaluated cell-wise on the diffuse band, with the source deposited as `m = V_Γ · w(F)` where

```
w(F) = F(1−F) / (√2·Cn)      ( = |∇F| on the equilibrium tanh profile )
```

Combined with the latent sink of G.4 (the *same* discrete field `m`), this is
enthalpy-conserving by construction and **self-pins** the band temperature to the liquidus: any
heat conducted into the band that the ice does not absorb is spent on melting, the interface
superheat scales as O(ε), and the flux-jump Stefan condition is recovered as ε → 0. The
solid-side conduction `κ_r,ice ∂θ/∂n|_sol` is accounted for automatically through the
phase-weighted diffusivity `κ_r(F)` (G.4) and the energy balance — `κ_r,ice` does **not** appear
in the melt law itself.

> **Why `w(F) = F(1−F)/(√2 Cn)` and not the measured `|∇F|`?** They coincide on the tanh band
> (same localization, same unit integral), but the CH relaxation is much slower than the melt
> deposition, so a `|∇F|`-weighted source keeps loading the nearly saturated tail cells (F → 1);
> the admissibility clip then silently destroys that melt mass *after* the temperature field paid
> its latent heat. In the A.3 gate this appeared as a ~10% enthalpy leak and a 4.4% slow front
> (λ = 0.325 instead of 0.340). `w(F)` vanishes quadratically at the tails, so the deposition
> itself never pushes cells out of [0,1].
>
> **Conservative admissibility clip (second leak channel).** Even with `w(F)`, the CH transient
> overshoots `C_L > 1` slightly behind the advancing front; the per-stage `[0,1]` clip removed
> ~5% of the influx. Under `PHASE_CHANGE` the clip is therefore made conservative in
> `VOF_DIFFUSE_step`: the clipped mass (returned by `diffuse_bound_liquid_fraction`) is
> redistributed back into the interfacial band by `diffuse_redistribute_liquid_mass_delta`
> (MPI-collective, F(1−F)-weighted, headroom-capped). With both elements the enthalpy budget
> closes and the front follows the Neumann law.

> **Why not the one-sided flux-jump evaluation** (`V_Γ = (St/Pe_T)([∂θ/∂n]_liq −
> κ_r,ice[∂θ/∂n]_sol)`, an earlier draft of this section)? On a diffuse band with *nothing pinning
> the interface temperature*, θ smooths across the band, the liquid- and solid-side one-sided
> gradients become equal, and their difference — hence the melt rate — collapses to ~0; the front
> stalls while heat leaks into the ice. This was observed directly in the A.3 test. The superheat
> form provides the pinning and the front speed in one energy-consistent term, so no separate
> liquidus-relaxation step is needed either.

- **`Λ*` (liquidus) ≠ `eos_Tmd_slope` (EOS density-maximum shift).** They look alike (both ∝ s)
  but are different physics: `Λ*` uses the seawater freezing slope `m_freeze`; `eos_Tmd_slope` uses
  the EOS coupling `cS`. Keep them separate.
- **Timestep note:** the melt source is explicit in the RK3 stages; its temperature-relaxation
  rate is `|∇F|/(Pe_T·ε) ~ 1/(4·Cn·Pe_T·ε)` at the band centre, so choose
  `dt ≲ Pe_T·ε·4Cn` (the A.3 test uses dt = 2e-4 with Pe_T = 10, ε = Cn ≈ 0.0029).
- **→ implemented by** A.1.4 (`stefan`, `T_melt`, `liquidus_slope`, `melt_band_eps`).

### G.4 Temperature with latent heat

```
∂θ/∂t + ∇·(u θ) = (1/Pe_T) ∇·( κ_r(F) ∇θ )  +  (1/St) (∂F/∂t)|_melt
```

- `Pe_T = Re·Pr = 10⁴` (reference); `κ_r(F) = κ(F)/κ_w = F + (1−F)κ_r,ice`: phase-weighted thermal
  diffusivity (heat diffuses in both phases).
- `(1/St)(∂F/∂t)|_melt`: latent sink — heat consumed where ice melts (the `V_Γ` source of G.2,
  **not** advective ∂F/∂t). `1/St = L/(cp ΔT) = St_Yang = 4`.
- **→ implemented by** A.1.3 (`CONC_VOF_PHASEWEIGHTED`, `kappa_ice_ratio`) for `κ_r(F)`; A.1.4
  (`Conc_add_source_RHS`) for the latent term.

### G.5 Salinity (salt confined to the liquid)

```
∂s/∂t + ∇·(u s) = (1/Pe_S) ∇·( F ∇s )
```

- `Pe_S = Re·Sc = 10⁶` (reference): the `Sc = 1000` salt field sets the 1440-cell resolution.
- Diffusivity masked by `F` ⇒ no salt flux into ice; meltwater dilution (`S = 0` released at the
  receding front) emerges from species conservation at the moving boundary.
- **→ implemented by** A.1.3 (`CONC_VOF_PHASEWEIGHTED`).

### G.6 Nonlinear equation of state (Roquet 2015, quadratic)

The linear `b = Σ Riᵢ cᵢ` is replaced by the nondimensional density anomaly

```
b(θ,s) = ρ'/Δρ_T = − β_T ( θ − T_md0 − slope·s )^q  +  β_S s ,     q = 2
```

which is the nondimensional form of `ρ' = −(Cb/2)(T − T0 − cS S)² + b0 S` with, for the reference
case, `β_T = 1`, `β_S = ΛT = 2 b0 Sm/(Cb ΔT²) = 1.75`, `T_md0 = (T0−Ti)/ΔT = 0.2`,
`slope = cS Sm/ΔT = −0.0625`. The salinity-shifted density maximum `θ_md(s) = T_md0 + slope·s`
captures the cold-water anomaly that drives the layering.

- **→ implemented by** A.1.2 (`EOS_NONLINEAR`; inputs `eos_q, eos_betaT, eos_betaS, eos_Tmd0,
  eos_Tmd_slope`).

### G.7 Dimensionless groups and how they enter the inputs

| Group | Definition | Reference value | Input |
|---|---|---|---|
| `Re` (free-fall) | `√(RaT/Pr)` | 1000 | `[flow] Re` |
| `Pe_T` | `Re·Pr` | 1e4 | `[conc] Pe[0]` |
| `Pe_S` | `Re·Sc` | 1e6 | `[conc] Pe[1]` |
| `St` | `cp ΔT/L` (= 1/St_Yang) | 0.25 | `[phase_change] stefan` |
| `β_S` | `ΛT = RaS/RaT` | 1.75 | `[eos] eos_betaS` |
| `Λ*` | `m_freeze Sm/ΔT` | ~0.014 | `[phase_change] liquidus_slope` |

### G.8 Method comparison with Yang / AFiD (added 2026-07-28)

**Provenance.** The PARTIES column is verified line-by-line against the source. The AFiD column is
**second-hand**: the arXiv letter (`docs/papers/yang2023_arxiv_source/main.tex`) contains no
phase-field equation at all — it defers the whole numerical method to Supplementary Material that
is not in the e-print bundle. Those entries come from the JFM text and the cited method papers
(Liu et al. 2021; Hester et al. 2020; Yang et al. 2022) and should be treated as such.

| | PARTIES (verified) | Yang / AFiD (second-hand) |
|---|---|---|
| Grids | **single grid**, everything at `N²` | **dual grid**: `u,T` on 288²; `S,φ` on 1440² (5×) |
| Phase eq. | conservative **Cahn–Hilliard**, `∂F/∂t + ∇·(uF) = Pe_CH⁻¹∇²ψ + m` | non-conservative **Allen–Cahn** |
| `ψ` | `W'(F) − Cn²∇²F`, `W' = ½F(1−F)(1−2F)` (`VOF_DIFFUSE.c:875`) | Gibbs–Thomson/superheat forcing, `εΔT/γ = 1` |
| Mobility | **constant**, `bih_coeff = Cn²/Pe_CH` (`VOF_DIFFUSE.c:1052`) | `D = 1.2κ` |
| Interface width | `Cn = 0.75Δx`; 5–95 % band `= 8.33 Cn` | `ε = Δx_refined = 1/1440`; band `= 5.89 ε` |
| Front motion | **added source** `m = (St/Pe_T ε)(θ−θ_L)·F(1−F)/(√2 Cn)` (`VOF_DIFFUSE.c:1299`) | intrinsic to the AC relaxation |
| Boundedness | needs `[0,1]` clip + conservative redistribution | naturally bounded by AC |
| Salt | `∂s/∂t + ∇·(us) = Pe_S⁻¹∇·(F∇s)`; conserves `∫s` | `Pe_S⁻¹(F+δ)⁻¹∇·[(F+δ)∇S] − S Ḟ/(F+δ)`; conserves `∫(F+δ)S` |
| Rigidity | implicit Brinkman diagonal `+ρ_f(2/τ)φ_s` (`msolve_cg.c:328`), `τ = 1e-4` | `η = dt` penalty **plus hard `u = 0` for `φ_s > 0.9`** |
| Time | semi-implicit low-storage RK3, `Γ={8/15,5/12,3/4}`, `Z={0,−17/60,−5/12}`, `B={4/15,1/15,1/6}` | AFiD RK3 |
| Space | MAC staggered, 2nd-order central; **WENO5** for `∇·(uF)`; HYPRE PCG+PFMG | 2nd-order staggered FD |

**The decisive number.** The 5–95 % interface band is `8.33 Cn` for PARTIES and `5.89 ε` for AFiD:

| | PARTIES 288² | PARTIES 720² | PARTIES 1440² | Yang |
|---|---:|---:|---:|---:|
| band | 2.17e-2 | 8.68e-3 | **4.34e-3** | **4.09e-3** |
| vs Yang | 5.3× | 2.1× | **1.06×** | — |

**At 1440² the two interface descriptions are the same width to 6 %.** The methods are therefore
*not* far apart on interface representation. What differs is that Yang feeds that sharp interface
from a temperature field resolved **5× more coarsely**, whereas PARTIES resolves `T` and `F` on the
same mesh. Since the `Sm=0` case has no salinity, their refinement buys nothing there and 288² sets
the heat delivery — which is exactly what §A.5.2 measured (Yang ≡ PARTIES at effective `N = 257`).

**Differences that are NOT drivers, with evidence:**

- **Brinkman vs direct forcing.** `δ = √(ντ) = 3.16e-4` is 0.46 cells at 1440² and 0.09 at 288² —
  sub-cell either way, so PARTIES' porous layer is already an effective no-slip wall. Consistent
  with the τ-sweep null (τ = 2.5e-4 closed 0.8 % of the gap).
- **Mobility model.** `Pe_CH × 5` (= Yang's `D` scale) moved `V` by 1.4e-4 at `t=40`.
- **Salt operator.** The Le=100 gate ran both against an exact binary-Stefan similarity solution:
  PARTIES' masked form was **1.5× better** on `λ` and 2–6× better on liquid-`S` `L∞`.

**Differences that DO matter:**

1. **Dual-grid vs single-grid** — explains the freshwater result (§A.5.2), essentially settled.
2. **CH vs AC boundedness.** CH conserves `F` but can be driven outside `[0,1]` by the source, so it
   *requires* the clip + redistribution; AC cannot. This is a genuine structural cost of our choice
   and it is what broke Step C (§A.5.5 item 3): off the tuned `Cn = 0.75Δx` mapping the clip
   destroys melt mass whose latent heat was already paid, leaking enthalpy at 1e-2.
3. **`κ` in the ice.** `kappa_ice_ratio[0] = 1.0` makes heat diffuse in ice at the *water* rate
   (physically `κ_ice/κ_water ≈ 8`). Low impact here because the ice is initialised at, and pinned
   to, the melting temperature — but it is an unverified modelling choice, not a match to Yang.

---

## STAGE A — YANG BENCHMARK

### A.1 Code modifications (all 3D)

Files referenced from the existing capability audit. Each item lists **file → function**, the work,
the **new flag**, and the **new inputs**.

#### A.1.1 Enable CONC + VOF_DIFFUSE together  ✅ already working on this branch
- The combination was validated by the pre-existing `testcases/TEST_VOFDIFFUSE_CONC_BOUSSINESQ`
  and `Canonical_VOFDIFFUSE_CONC_BOUSSINESQ` cases (CONC + VOF_DIFFUSE + BOUSSINESQ, no
  particles): allocation, ghost exchange, and the one-fluid velocity for `Conc` advection all
  coexist without changes to `Cart3d_initialize_primitive_data`.
- **Initialization-order caveat** (matters for A.1.6): `Cart3d_initialize_primitive_data`
  initializes the `Conc` fields **before** the VOF field, so scalar init types must not read `F` —
  the new inits 30/31 recompute the slab profile analytically from `vof_slab_x0`/`Cn`.
- **RK3 substep order (verified in `Temporal_int_all_the_equations`)**: `VOF_DIFFUSE_step`
  (melt rate from `F^{k-1}`, `θ^{k-1}` computed at its top → CH source) → `Conc_int_equations`
  (latent sink, same stage, same RK weights) → `Velocity_add_buoyancy_2_RHS` (A.1.2).
- **No new flag** (uses existing `CONC`,`VOF_DIFFUSE`). **`NConc = 2`** (field 0 = T, field 1 = S).

#### A.1.2 Nonlinear EOS (Roquet 2015, quadratic q=2)  ✅ implemented and validated 2026-07-05
- **`PARTIES/src/Eulerian/Velocity.c`** (`Velocity_add_buoyancy_2_RHS`) — replaces the
  linear `b = Σ Riᵢ cᵢ` with
  `b(θ,s) = −β_T·|θ − T_md0 − slope·s|^q + β_S·s` (Boussinesq, ρ₂/ρ₁ = 1 for ice/water; the
  anomaly is symmetric about the density maximum, so the exponent acts on the absolute value —
  identical to `(·)²` for q = 2). Requires `NConc ≥ 2` (field 0 = θ, field 1 = s; runtime abort
  otherwise).
- **Inputs are parsed via `src/IO/default.inp`** (the X-macro registry used by `Input.c`) — new
  `[eos]` block keys there, defaults inert.
- **New flag:** `EOS_NONLINEAR` (requires `BOUSSINESQ` + `CONC`; enforced in `Boundary.h`).
- **New inputs (`[eos]`):** `eos_q`, `eos_betaT`, `eos_betaS`, `eos_Tmd0`, `eos_Tmd_slope`.

#### A.1.3 Phase-weighted scalar diffusivity (heat in both phases; salt only in liquid)  ✅ implemented
- **`PARTIES/src/Eulerian/Conc.c`** — both halves of the Crank–Nicolson split get the same
  face-averaged coefficient `κ(F)/κ_liq = F + (1−F)·kr` (helper `Conc_phase_kappa`):
  `Conc_set_conv_viscous_central_mixed` (explicit half) and `Conc_laplacian` (implicit CG
  operator, now variable-coefficient). Temperature uses a finite `kr`; salinity uses `kr = 0`
  (salt flux vanishes in ice, so S stays 0 inside ice with no melting).
- **New flag:** `CONC_VOF_PHASEWEIGHTED` (requires `VOF_DIFFUSE` + `CONC` + `CONC_CENTRAL` +
  `CONC_FULLY_IMPLICIT`; incompatible with `VAR_VISC`/`VOF_SCALAR`; enforced in `Boundary.h`).
- **New input (`[conc]`):** `kappa_ice_ratio` — a **per-field array** like `Pe`, e.g.
  `kappa_ice_ratio = {1.0, 0.0}` for Yang (T conducts in ice, S masked). Default `{1, 0}`.

#### A.1.4 Phase change — Stefan source in CH + latent heat in T  ✅ implemented
- **`PARTIES/src/Eulerian/VOF_DIFFUSE.c`** — new function `VOF_DIFFUSE_compute_melt_rate`,
  called at the top of `VOF_DIFFUSE_step` every RK stage (from `F^{k-1}`, `θ^{k-1}`): superheat
  melt rate `V_Γ = (St/Pe_T)·(θ − θ_L(s))/melt_band_eps` (see G.3 — the one-sided flux-jump form
  of the earlier draft stalls on an unpinned diffuse band and was replaced); localized source
  `m = V_Γ·F(1−F)/(√2·Cn)` added to the explicit CH RHS (the tail-safe equivalent of `V_Γ|∇F|`,
  see G.3), so it rides the same GAMMA/ZETA RK3 combination and `ch_rhs_n` history as
  advection/diffusion; `m` is kept in `vof->melt_src` (+ `melt_src_old` for the ZETA history term
  of the temperature equation). **[MPI]** the melt law is fully pointwise (no stencil), so it
  needs no halo refresh and is rank-count independent by construction; `θ`/`s` halos are refreshed
  at the top of `Temporal_int_all_the_equations` anyway. Ice-volume `V(t)` style diagnostics must
  use `MPI_Allreduce`.
- **`PARTIES/src/Eulerian/Conc.c`** — new `Conc_add_latent_heat_RHS` (called from
  `Conc_int_equations` for field 0, after `Conc_set_RHS`): latent sink
  `−St⁻¹·(GAMB_k·m + ZETB_k·m_old)`, i.e. the *same* discrete melt field with the *same* RK
  weights as the CH source — the discrete enthalpy `∫θ + St⁻¹∫F` then changes only through
  boundary fluxes.
- **`VOF_DIFFUSE_step`** — under `PHASE_CHANGE` the per-stage `[0,1]` admissibility clip is made
  **conservative**: the clipped mass is redistributed back into the interfacial band (see the
  G.3 note; otherwise the clip is a hidden enthalpy leak that slows the front).
- **Liquidus:** enters directly through `θ_L(s) = T_melt − liquidus_slope·s` in the melt law; the
  negative feedback pins the band temperature, so **no separate pinning/relaxation step exists**.
  **Note:** `liquidus_slope` is the **freezing-point depression**, a *different* quantity from the
  EOS coupling `eos_Tmd_slope` — do not reuse `cS` here.
- **New flag:** `PHASE_CHANGE` (requires `VOF_DIFFUSE` + `CONC`; enforced in `Boundary.h`).
- **New inputs (`[phase_change]`):** `stefan` (= cp·ΔT/L; 0 disables melting), `T_melt`,
  `liquidus_slope`, `melt_band_eps` (melt-law band scale ε; ≤ 0 → defaults to `Cn`).

#### A.1.5 Ice rigidity — volume penalization (Darcy damping)  ✅ implemented and validated 2026-07-05
- **`PARTIES/src/Eulerian/lsolver/msolve_cg.c`** (`matVec`, `vel_jacobi_diag`) — the Darcy
  damping is added as a **positive local diagonal** in the implicit velocity operator:
  `+ 2·ρ_f·φ_s,f/darcy_tau` (the factor 2 matches the code's CN convention where explicit sources
  carry a 2×). Face-centered ice mask `φ_s,f = 1 − F_f` (F = `C_L` under `VOF_DIFFUSE`), clipped
  to `[0,1]`. The operator stays SPD; the Jacobi preconditioner gets the same diagonal. This is
  the production path; it avoids the explicit stability limit `Δt_stage ≲ τ_p`.
- **`PARTIES/src/Eulerian/lsolver/msolve_direct.c`** (`Velocity_solve_explicit`) — the same
  diagonal denominator in the explicit intermediate velocity used to estimate IBM forcing:
  `u_hat = RHS / (ρ_f/(α_kΔt) + 2ρ_fφ_s,f/darcy_tau)` (active whenever `ICE_PENALIZATION` is on,
  which covers the later `LAG_PARTICLE_RESOLVED` Stage-B use).
- **`PARTIES/src/Eulerian/Velocity.c`** — do **not** add the production term as a plain explicit
  RHS force. A weak explicit version is acceptable only as a debug switch for `darcy_tau` no smaller
  than the RK substep.
- A large `mu2/mu1` can be used as a debugging rigidity backup, but the Yang production input uses
  `mu2/mu1 = 1` because the paper's momentum equation has constant viscosity; production rigidity
  comes from `ICE_PENALIZATION`. A useful first `darcy_tau` estimate is
  `darcy_tau ≈ Re·Δx²`, which gives a
  Brinkman damping layer `sqrt((1/Re)·darcy_tau)` of order one grid cell. Then run a sensitivity
  sweep and choose the largest damping time that makes ice velocity leakage irrelevant.
- **New flag:** `ICE_PENALIZATION`.
- **New input (`[vof]` or `[phase_change]`):** `darcy_tau`.

#### A.1.6 New initial condition — vertical ice slab + T step + S linear profile  ✅ implemented
- **`PARTIES/src/IO/VoF_Init.c`** (dispatched from `Cart3d.c`) — new VOF `init_type = 27`
  (`VoF_init_vertical_ice_slab`): `F = 0` (ice) for `x ≥ vof_slab_x0·Lx`, `F = 1` (water)
  elsewhere, tanh-smoothed over the CH band (`F = 0.5(1 − tanh((x − x_int)/(2√2·Cn)))`). The
  interface position is a new input `vof_slab_x0` (**default 0.9** = Yang; the A.3 Stefan test
  sets it to 0.05 for a thin water film at the hot wall — additive, so Yang is unchanged).
- **New input (`[vof]`):** `vof_slab_x0` (interface location as a fraction of `Lx`; default 0.9).
- **`PARTIES/src/IO/Initial_Conditions.c`** — new `conc_init_type`:
  - `30` (T, `Conc_init_ice_slab_T`): `θ = cbd0` (default 1) in water, `θ = theta_ice`
    (default 0; set < 0 for the subcooled two-phase Neumann variant) in ice.
  - `31` (S, `Conc_init_ice_slab_S_ylinear`): `s = s_bot + (s_top − s_bot)·y/Ly` in water,
    `s = 0` in ice, with `s_top = cbd2`, `s_bot = cbd5`.
  Both blend with the **same analytic tanh slab profile** as init 27 (computed from
  `vof_slab_x0`/`Cn`, *not* by sampling F — `Cart3d_initialize_primitive_data` initializes the
  Conc fields **before** the VOF field, so the scalar inits must not read F).
- **New inputs (`[conc]`):** `theta_ice` (ice temperature for init 30); reuse `cbd0` = θ_water,
  `cbd2` = s_top, `cbd5` = s_bot.

> **Checklist (Stage A code) — all implemented:** `EOS_NONLINEAR`, `CONC_VOF_PHASEWEIGHTED`,
> `PHASE_CHANGE`, `ICE_PENALIZATION`, init types 27 / 30 / 31 (init 27 takes `vof_slab_x0`,
> default 0.9, for the A.3 Stefan gate), RK3 substep ordering, `[eos]` + `[phase_change]` input
> parsing. No IBM, no particles, no surface tension.

### A.2 Debugging pass (tiny reduced domain — code shake-out only, NOT validation)  ✅ PASSED 2026-07-05

> **Result:** the reduced integrated Yang checkout (`PARTIES/testcases/ECCO_TESTS/Yang_checkout/`,
> `analyze_yang_checkout.py`, SLURM job 273214) ran the **full Yang flag set**
> (`BOUSSINESQ + EOS_NONLINEAR + CONC_VOF_PHASEWEIGHTED + PHASE_CHANGE + ICE_PENALIZATION`,
> `NConc = 2`) on 128×128×4 (1 and 4 ranks), 128×128×8, and a static control — 50 steps each,
> all 11 gate checks PASS: no NaN/Inf, salt drift 3.7·10⁻¹⁶, enthalpy drift 1.5·10⁻⁸,
> 1-vs-4-rank worst field difference 1.3·10⁻¹³, NZ=4 vs NZ=8 z-mean difference 2.0·10⁻⁹,
> static control exactly quiescent with front displacement 0.0012 cells. The first launch
> exposed (and fixed) an off-by-one pressure-loop abort in `Temporal_int.c` — the projection
> now aborts only when the residual is still above tolerance after the final iteration.
> Full record: `PARTIES/testcases/ECCO_TESTS/Yang_checkout/YANG_CHECKOUT_REPORT.md`.

Purpose: confirm it **compiles, runs, exchanges ghosts/MPI in 3D, and conserves** — nothing
physical is being validated. Run a handful of these, each only tens of steps.

> **Sequencing as-built:** all four quantitative gates A.3–A.3d passed independently first.
> They isolate individual modules and do not replace this integrated shakeout, which has now
> also passed — **the A.4 production run is unblocked**.

```ini
# DEBUG grid: 10x coarser, reduced Sc so salt is representable on a coarse mesh.
[geometry]
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0
zmin = 0.0
zmax = 0.03125          # 4 cubic cells

[grid]
NXM = 128
NYM = 128
NZM = 4

[simulation]
time_max = 0.05         # 50 fixed steps at dt=0.001; extend only if diagnostics need it
default_dt = 0.001
output_time_interval = 0.5

[flow]
Re = 1000.0

[conc]
NConc = 2
Pe = {1e4, 1e4}         # Sc=10 (Pe_S=Re*10) for cheap debugging only
# ...all Stage-A flags/keys as in A.4, with the values below.
```

**What to check (diagnostics, not physics):** (1) no NaN/Inf over 50 steps; (2) S stays exactly 0
inside the ice (no salt leak through the CH interface); (3) global heat + salt budgets close to
round-off when melting is toggled off; (4) interface stays put with `darcy_tau` large and melting
off; (5) implicit CG (viscous + CH) converges at the chosen `mu2/mu1`; (6) one- and four-rank
runs agree within solver tolerance; (7) the periodic-direction mean buoyancy is independent of
the thin-grid count (`NZM=4` versus 8), guarding the A.3c/A.3d endpoint fix. Fix bugs here; do
**not** interpret the coarse/low-Sc result as a benchmark.

### A.3 Preliminary validation — 1D Stefan (Neumann) problem  ✅ PASSED 2026-07-02

> **Result:** λ = 0.3400 vs exact 0.3401 (0.03 % error), similarity-profile error ≤ 2.5·10⁻³,
> exact quiescence, enthalpy closure 6.9·10⁻⁵, machine-precision 1-vs-4-rank consistency.
> Full record: [STAGE_A_STEFAN_GATE_REPORT.md](STAGE_A_STEFAN_GATE_REPORT.md); working case:
> `PARTIES/testcases/1DStefan/` (`analyze_stefan.py`, exit 0 = pass).

**Purpose.** The first *quantitative* physics check. It
certifies the new melting kernel — the Stefan source in the CH equation (G.2), the Stefan
condition / interface pinning (G.3), and the latent-heat sink in `θ` (G.4) — against
the **exact Neumann analytical solution**, in isolation. It is cheap (1D thin strip, low `Pe_T`,
seconds–minutes) and de-risks the 8.3 M-cell Yang run: a Yang discrepancy cannot tell you whether
the bug is in the Stefan kernel, the EOS, the salt coupling, or the penalization, whereas this test
pins the kernel down alone. It is also the place where the **`stefan = cpΔT/L` vs `St_Yang =
L/cpΔT` inversion** (see G.7) is caught — a wrong sense makes the front advance at a visibly wrong
√t rate.

**Validates which code modules:** A.1.4 (`PHASE_CHANGE`: melt rate `V_Γ`, the `V_Γ|∇F|` CH
source, the `1/St` latent sink) and A.1.3 (`CONC_VOF_PHASEWEIGHTED`: `κ_r(F)` — exercised
fully by the optional two-phase variant). It deliberately does **not** exercise A.1.2 (EOS), A.1.5
(penalization), or the salt field — those are off.

#### The problem and its analytical solution

One-phase melting in a semi-infinite domain: hot wall at `x = 0` held at `θ_w = 1`; ice initially at
the melt temperature `θ_m = 0` (so the solid carries no gradient — one-phase). The liquid conducts,
the front `X(t)` recedes into the ice. With **buoyancy off the velocity is identically zero** (see
the flags below and the velocity argument under G.1), so G.4 collapses to pure conduction
plus the latent sink, and the nondimensional thermal diffusivity is `α = 1/Pe_T`:

```
θ(x,t) = 1 − erf( x / (2√(α t)) ) / erf(λ) ,     α = 1/Pe_T

X(t)   = 2 λ √( t / Pe_T )

λ e^{λ²} erf(λ) = St / √π ,     St = stefan · (θ_w − θ_m) = stefan
```

For the values below (`stefan = 0.25`, `θ_w − θ_m = 1`): `St = 0.25` ⇒ `St/√π = 0.141` ⇒
**`λ ≈ 0.340`**. The front therefore follows `X(t) = 0.680 √(t/Pe_T)`.

**Optional two-phase variant** (run *after* the one-phase gate passes, to also certify the
heat-in-ice term `κ_r,ice`): subcool the ice to `θ_∞ < 0`. With equal properties
(`kappa_ice_ratio = 1.0`, `α_s = α_l`) the front constant solves

```
St_l / (e^{λ²} erf(λ))  −  St_s / (e^{λ²} erfc(λ))  =  λ √π ,
St_l = stefan·(θ_w − θ_m),   St_s = stefan·(θ_m − θ_∞)
```

This is the only place the solid-side `∂θ/∂n` branch of `V_Γ` and the `κ_r(F)` ice conduction are
checked against an exact answer. (One-phase uses conc_init_type 30 (A.1.6) as-is; the two-phase
variant needs the ice value `θ_∞` exposed as an input — a small additive generalization of init 30.)

#### `Boundary.h` (Stefan test)

3D (thin strip exercises the 3D path cheaply). `x` = no-slip walls (Dirichlet-`θ` hot wall at
west, see BCs); `z` = **periodic** and thin; `y` = **no-slip walls with no-flux `θ`** (identical
to periodic for this y-uniform 1D problem, and it stays on the well-tested wall code path —
`YPERIODIC` + scalars is untested). Define `VOF_DIFFUSE`, `CONC`, `PHASE_CHANGE`,
`CONC_VOF_PHASEWEIGHTED`, `USE_HYPRE`. **Leave undefined** `EOS_NONLINEAR`, `BOUSSINESQ`,
`ICE_PENALIZATION`, `SURFACE_TENSION`, `VOF_GRAVITY`, `VOF_IBM`, `LAG_PARTICLE_RESOLVED`,
`PARTICLE_RELEASE`, `TWOD_CARTESIAN`. `NConc = 1` (temperature only — **no salt**).

> **Why buoyancy must be off.** The melt front is a *vertical* interface, so `∇θ` is horizontal; a
> horizontal density gradient under gravity has no hydrostatic balance and would drive baroclinic
> convection that destroys the 1D similarity solution. With `EOS_NONLINEAR`/`BOUSSINESQ` undefined
> the momentum equation has no body force; the remaining viscous (and, if ever enabled, Darcy)
> terms are purely dissipative, and ice/water share `ρ₀` so melting causes no expansion flow ⇒ the
> velocity stays at round-off. `max|u| < 10⁻¹⁰` is therefore itself a pass check (a rise signals
> leaked buoyancy or a momentum-injecting melt source).

#### Inputs (`parties.inp` — the PARTIES input file name; full working case in
`PARTIES/testcases/1DStefan/`)

```ini
[geometry]
# thin strip: front moves in x; y,z thin (z periodic, y no-flux walls)
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 0.015625         # 4 cells
zmin = 0.0
zmax = 0.015625         # 4 cells

[grid]
NXM = 256               # dx = 1/256 ; ~3-4 cells across the CH band
NYM = 4
NZM = 4

[simulation]
# X(t)=0.680*sqrt(t/Pe_T); with Pe_T=10 the front reaches x~0.8 by t~14
time_max = 15.0
default_dt = 2.0e-4     # explicit melt-source relaxation limit: dt < Pe_T*eps*4Cn ~ 5e-4
max_dt = 2.0e-4
constant_dt = 1         # u = 0, so the advective CFL is meaningless here
cfl = 0.3
output_time_interval = 0.5

[flow]
Re = 1.0                # irrelevant (u=0); keeps viscous CG well-scaled
vel_init_type = 0       # quiescent

[vof]
rho1 = 1.0              # water = phase1 (F=1)
rho2 = 1.0              # ice   = phase2 (F=0); Boussinesq rho2=rho1
mu1 = 1.0
mu2 = 100.0            # rigidity backup only (no flow anyway)
init_type = 27         # vertical ice slab (reads the interface position)
vof_slab_x0 = 0.05     # thin water film at the hot wall; ice for x >= 0.05*Lx (default 0.9 = Yang)
Cn = 0.0029296875      # 0.75*dx ; CH band ~3-4 cells
# darcy_tau / ICE_PENALIZATION OFF

[phase_change]
stefan = 0.25          # = cp*dT/L ; St = stefan*(theta_w-theta_m) = 0.25 -> lambda~0.340
T_melt = 0.0           # theta of melting at S=0
liquidus_slope = 0.0   # no salt -> no freezing-point depression
melt_band_eps = 0.0029296875   # melt-law band scale (= Cn)

[conc]
# field 0 = T only
NConc = 1
Pe = {10}              # Pe_T small so the front moves fast and cheap (no salt Batchelor limit)
kappa_ice_ratio = {1.0} # per-field array: heat conducts in ice (exercised by the two-phase variant)
conc_init_type = {30}  # theta = cbd0 in water, theta_ice in ice  (one-phase IC)
cbd0 = 1.0             # theta_w (water)
theta_ice = 0.0        # theta_inf (ice); set < 0 for the two-phase variant
# Scalar BC convention in the code: A*dtheta/dn + B*theta = C on each face
# (A multiplies the GRADIENT — an earlier draft of this block had A and B swapped).
# Hot wall: Dirichlet theta = 1 on the WEST (x=0) face:
BC_AW = {0}
BC_BW = {1}
BC_CW = {1}
# All other non-periodic faces: no-flux (A=1, B=0, C=0); east stays cold (semi-infinite)
BC_AE = {1}
BC_BE = {0}
BC_CE = {0}
BC_AN = {1}
BC_BN = {0}
BC_CN = {0}
BC_AS = {1}
BC_BS = {0}
BC_CS = {0}
```

> For the **two-phase variant**, set the ice far-field value `theta_ice < 0` (e.g. `−0.5`), keep
> the east wall no-flux, and compare against the two-phase `λ` above.

#### What to check (the actual gate)

1. **Front law.** Extract the `F = 0.5` location `X(t)` and confirm it collapses onto
   `X(t) = 2λ√(t/Pe_T)` with `λ ≈ 0.340` (use a virtual time origin `t₀ = (X₀/0.680)²·Pe_T` from the
   initial film `X₀ = 0.05·Lx`, i.e. compare for `t > t₀`).
2. **Similarity profile.** `θ(x,t)` collapses onto `1 − erf(η)/erf(λ)` when plotted vs
   `η = x/(2√(αt))` at several times.
3. **Diffuse → sharp convergence.** Refine `Cn`/`melt_band_eps` (and the grid); the measured `λ`
   converges toward the analytic value — this validates the regularized `V_Γ|∇F|` source and fixes
   defensible `Cn`/`melt_band_eps` for Yang (A.4).
4. **Convention.** A √t front with the *correct* `λ` confirms `stefan`, the `1/St` latent factor,
   and their signs are wired correctly (the inversion footgun of G.7).
5. **Quiescence.** `max|u| < 10⁻¹⁰` throughout (independent check that buoyancy is truly off and the
   melt source injects no spurious momentum).
6. **Enthalpy budget.** Total enthalpy `∫θ + (1/St)∫F` (sensible plus latent; equivalently
   `∫θ − (1/St)∫(1−F)` + const — note the sign: converting ice to water *stores* latent heat)
   balances the integrated hot-wall heat flux — the independent check on the latent coupling.
   Offline (from the saved fields) the closure is limited by the time-quadrature of the sparse
   outputs, not round-off.

**[MPI]** the `F = 0.5` front extraction and the enthalpy/heat-flux budgets use `MPI_Allreduce`, not
per-rank values; refresh the `F`/`θ` halos before the `∂θ/∂n` stencils exactly as in A.1.4.

**Pass criterion:** criteria 1–6 hold and `λ` converges to `0.340 ± O(Δx, Cn)`. Only then proceed
to A.3b. (As-built note: this gate ran *before* the A.2 shakeout — it needs only the Stefan-gate
flag subset (`CONC`, `VOF_DIFFUSE`, `PHASE_CHANGE`, `CONC_VOF_PHASEWEIGHTED`), not the full Yang
set; the only code change beyond A.1 was the `vof_slab_x0` input on init 27, default 0.9 so Yang
is unchanged.)

### A.3b Intermediate validation — 1D binary (salty) Stefan problem  ✅ PASSED 2026-07-03

> **Result (first attempt, no code changes needed):** λ = 0.25320 vs exact 0.25554 (**0.92 %
> error**, front-fit rms 4.0·10⁻⁵ ≈ Δx/49) — 23× closer to the correct value than to either
> broken-coupling reference; two-sided θ collapse ≤ 1.7·10⁻², s collapse ≤ 4.5·10⁻³;
> s(X) → 0.1666 and θ(X) → −0.0780 converging monotonically toward the exact 0.1722/−0.0861
> (O(band) sampling offset, shrinking as the layers widen); salt conserved to **2.9·10⁻¹³**,
> enthalpy to **3.0·10⁻¹⁰**, salt-in-ice ≤ 3.6·10⁻¹¹, velocity **exactly 0.0**.
> Rank-count independence re-confirmed for the coupled case (1 vs 4 ranks at t = 0.45:
> Δθ ≤ 2.1·10⁻¹⁵, Δs ≤ 1.0·10⁻¹⁵, ΔF ≤ 3.2·10⁻¹⁵, Δu = 0).
> Full record: [STAGE_A_SALTY_STEFAN_GATE_REPORT.md](STAGE_A_SALTY_STEFAN_GATE_REPORT.md);
> gate + plots: `PARTIES/testcases/1DStefanSalty/analyze_salty_stefan.py` (exit 0 = pass);
> results in `stefan_salty_results.csv` / `stefan_salty_timeseries.csv` /
> `stefan_salty_profiles.csv` and `fig_*.png` in the same folder.

**Purpose.** The A.3 gate certified the melting kernel for pure water; the Yang run additionally
depends on the **salt–melt coupling**, which A.3 deliberately left off. This gate certifies that
coupling against an exact analytical solution, still in a cheap quiescent 1D strip:

1. the **salinity-dependent liquidus** in the melt law, `θ_L(s) = T_melt − Λ*·s` (A.1.4) — with
   the depression made an **O(1) player** (`liquidus_slope = 0.5` instead of Yang's 0.014) so a
   miswired coupling shifts the front constant by 15 %, not 1 %;
2. the **F-masked salt diffusivity** `kr = 0` (A.1.3): no salt flux into ice, and the **meltwater
   dilution boundary condition emerging from species conservation** at the receding front — the
   central claim of G.5, never before checked quantitatively;
3. the **solid-side conduction path** `κ_r(F)` with `kappa_ice_ratio = 1` (A.1.3): with the
   interface depressed to `θ_i < 0` and the ice at `θ = 0`, the solid conducts heat *into* the
   front — the two-phase Neumann branch (the "optional variant" of A.3, exercised here with the
   sign it actually has in salty melting);
4. the **`NConc = 2` field indexing** (0 = θ, 1 = s) and per-field `Pe`/`kappa_ice_ratio`/BC
   arrays — exactly the Yang configuration.

It deliberately does **not** exercise A.1.2 (`EOS_NONLINEAR`/buoyancy) or A.1.5
(`ICE_PENALIZATION`) — a quiescent 1D similarity solution requires zero flow, so those have no
exact reference and are covered by the A.2 shakeout diagnostics and the Yang comparison itself.

#### The problem and its analytical solution

Fresh ice and warm salty water, both effectively semi-infinite, no walls driving anything:
liquid (`θ = θ_∞ = 1`, `s = s_∞ = 1`) for `x < x_int`, ice (`θ = θ_ice = 0`, `s = 0`) for
`x > x_int`. All scalar BCs no-flux, `u = 0` identically (no buoyancy, equal densities). Melting
is driven by the initial liquid superheat; the front `X(t)` advances into the ice as `√t`. Because
uniform initial states in each phase *are* the `t → 0` limit of the similarity solution, the
similarity law holds essentially from the start (only the diffuse-band relaxation introduces a
small virtual time origin).

With `η = (x − x_int)/(2√(αt))`, `α = 1/Pe_T`, and `η_s = (x − x_int)/(2√(Dt))`, `D = 1/Pe_S`:

```
X(t) = x_int + 2λ√(αt)

liquid (η < λ):   θ = θ_∞  + (θ_i − θ_∞ )·(1 + erf(η)) /(1 + erf(λ))
solid  (η > λ):   θ = θ_ice + (θ_i − θ_ice)·erfc(η)/erfc(λ)
liquid (η_s < λ_s): s = s_∞ + (s_i − s_∞)·(1 + erf(η_s))/(1 + erf(λ_s)) ,   λ_s = λ√Le
```

The three interface conditions close the system for `(λ, s_i, θ_i)`:

```
liquidus:        θ_i = T_melt − Λ*·s_i

salt dilution:   (1/Pe_S)·∂s/∂x|_i = −Ẋ·s_i   (fresh melt dilutes the front)
                 ⇒  s_i = s_∞ / [ 1 + √π·λ_s·e^{λ_s²}·(1 + erf(λ_s)) ]

Stefan:          (1/St)·Ẋ = (1/Pe_T)·( ∂θ/∂x|_solid − ∂θ/∂x|_liquid )
                 ⇒  λ√π e^{λ²}/stefan = (θ_∞ − θ_i)/(1 + erf(λ)) + (θ_ice − θ_i)/erfc(λ)
```

(The A.3 two-phase relation is the wall-anchored special case; here both terms *add* because the
liquidus depression puts the interface below both far fields — freezing-point depression makes ice
melt **faster**, the road-salt effect.)

**Reference constants** for the inputs below (`stefan = 0.5`, `Λ* = 0.5`, `Pe_T = 100`,
`Pe_S = 1000`, `Le = 10`), from the coupled transcendental system (solved in
`analyze_salty_stefan.py`):

| Quantity | Value | Broken-coupling references |
|---|---|---|
| `λ` (front constant) | **0.25554** | 0.23391 (−8.5 %) if solid conduction lost; 0.21688 (−15.1 %) if liquidus/salt coupling lost |
| `λ_s = λ√Le` | 0.80810 | — |
| `s_i` (interface salinity) | **0.17220** | 1.0 if dilution BC lost |
| `θ_i = −Λ*·s_i` (interface temperature) | **−0.08610** | 0 if liquidus lost |

A ~1 % front-law measurement (A.3 achieved 0.03 %) separates the correct `λ` from either broken
path by > 8×, so the gate is sharply discriminating.

#### `Boundary.h` (salty Stefan gate) — identical to the A.3 gate

Define `VOF_DIFFUSE`, `CONC`, `PHASE_CHANGE`, `CONC_VOF_PHASEWEIGHTED`, `USE_HYPRE`;
**leave undefined** `EOS_NONLINEAR`, `BOUSSINESQ`, `ICE_PENALIZATION`, `SURFACE_TENSION`,
`VOF_GRAVITY`, `VOF_IBM`, `LAG_PARTICLE_RESOLVED`, `PARTICLE_RELEASE`, `TWOD_CARTESIAN`.
`x` = no-slip walls, `y` = no-slip walls, `z` = periodic (all as in A.3 — no
`Boundary.h` change is needed between the two gates; `NConc` is a runtime input).

#### Inputs (`parties.inp`; full working case in `PARTIES/testcases/1DStefanSalty/`)

```ini
[geometry]
# both phases must stay effectively semi-infinite: at t = 4.5 the liquid gap
# is 4.1 and the (shrinking) solid gap 3.0 diffusion lengths sqrt(alpha*t)
xmin = 0.0
xmax = 1.5
ymin = 0.0
ymax = 0.0078125        # 4 cells
zmin = 0.0
zmax = 0.0078125        # 4 cells

[grid]
NXM = 768               # dx = 1/512; CH band ~3-4 cells; salt BL 23 cells at t=0.5
NYM = 4
NZM = 4

[simulation]
time_max = 4.5          # front travels ~0.11 = 55 cells; far fields stay clean
output_time_interval = 0.15
default_dt = 2.0e-4     # melt-source limit Pe_T*eps*4Cn ~ 8.6e-4; CH explicit ~5e-4
max_dt = 2.0e-4
constant_dt = 1         # u = 0, advective CFL meaningless

[flow]
Re = 1.0
vel_init_type = 0       # quiescent

[vof]
rho1 = 1.0
rho2 = 1.0              # ice/water share rho0 (Boussinesq-style)
mu1 = 1.0
mu2 = 100.0
init_type = 27
vof_slab_x0 = 0.5       # interface at x = 0.75 (fraction of Lx = 1.5)
Cn = 0.00146484375      # 0.75*dx
Pe_CH = 614.4           # = 0.9/Cn, same convention as A.3

[phase_change]
stefan = 0.5            # St_l = stefan*(theta_inf - theta_i)
T_melt = 0.0
liquidus_slope = 0.5    # EXAGGERATED (Yang: 0.014) to make the coupling O(1)
melt_band_eps = 0.00146484375

[conc]
NConc = 2               # field 0 = theta, field 1 = s  (Yang layout)
Pe = {100.0, 1000.0}    # Pe_T, Pe_S  (Le = 10)
kappa_ice_ratio = {1.0, 0.0}   # heat conducts in ice; salt masked - Yang values
conc_init_type = {30, 31}
cbd0 = 1.0              # theta_infinity (water)
theta_ice = 0.0         # theta far field in ice
cbd2 = 1.0              # s_top = s_bot = 1 -> UNIFORM s in water (init 31)
cbd5 = 1.0
# ALL faces no-flux for BOTH fields (A=1, B=0, C=0): melting is driven by the
# initial superheat, so total enthalpy AND total salt are exactly conserved -
# two machine-precision budget checks with no wall-flux quadrature caveat.
BC_AW = {1.0, 1.0}
BC_BW = {0.0, 0.0}
BC_CW = {0.0, 0.0}
# ... identical no-flux triplets on E/N/S faces; z periodic.
```

#### What to check (the actual gate)

1. **Front law (primary).** `F = 0.5` location fits `X(t) = x_int + 2λ√(α(t − t₀))` over
   `t ∈ [0.5, 4.5]` with `λ = 0.2555 ± 3 %` — i.e. > 8σ away from both broken-coupling values.
2. **Temperature similarity (two-sided).** `θ(η)` collapses on the analytic profile in *both*
   phases (max error ≤ 0.05 outside the diffuse band) — certifies the solid-conduction branch.
3. **Salinity similarity.** `s(η_s)` collapses on the analytic dilution profile
   (max error ≤ 0.05·s_∞ outside the band).
4. **Interface values.** `s` and `θ` sampled at `X(t)` are time-independent and match
   `s_i = 0.1722` (±15 %) and `θ_i = −0.0861` (±0.03 absolute) — the tolerances are the
   half-band sampling uncertainty, O(band·|∇s|), not a physics allowance.
5. **Salt conservation.** `∫s dV` constant to ≤ 10⁻⁸ relative (no-flux everywhere, u = 0, source-
   free — any drift is a conservation bug in the F-masked operator).
6. **No salt in ice.** `max s` beyond the band (`x > X + 20Δx`) stays ≤ 10⁻⁸.
7. **Enthalpy conservation.** `∫(θ + F/St) dV` constant to ≤ 10⁻⁸ relative (all-no-flux ⇒ *total*
   conservation; strictly stronger than the A.3 wall-flux closure).
8. **Quiescence.** `max|u| < 10⁻¹⁰` throughout.

**Pass criterion:** all eight hold. The automated gate is
`PARTIES/testcases/1DStefanSalty/analyze_salty_stefan.py` (exit 0 = pass); it writes
`stefan_salty_results.csv` (time series + fit summary) and the `fig_*.png` plots. After this
gate, complete A.3c, A.3d, and the A.2 full-flag shakeout before A.4 (all done as of
2026-07-05).

**Resolution/parameter caveats.** The salt boundary layer `2√(Dt)` must stay well wider than the
CH band for the sharp-interface comparison to hold (here ≥ 6 bands at the fit-window start,
growing to ~17); the measured `s(X)`/`θ(X)` carry O(half-band × interface gradient) sampling
error (~10 % of `s_i` at `t = 1`), which is why criterion 4's tolerances are looser than
criterion 1's. `Le = 10` is a deliberate compromise: large enough that the salt layer is
distinctly thinner than the thermal layer (`√Le ≈ 3.2×`), small enough to resolve on a 1D strip —
the Yang-scale `Le = 100` steepening is a resolution question already handled by the 1440² grid
choice (A.4), not a coupling question. As in A.3, refine `Cn`/`melt_band_eps`/grid together if
the measured `λ` needs to be pushed closer to the analytic value.

### A.3c Intermediate validation — nonlinear EOS reconstruction  ✅ PASSED 2026-07-05

**Purpose.** Certify A.1.2 independently of phase change and ice drag. The all-liquid test uses
frozen, one-dimensional `(theta,s)` profiles and gravity along the thin periodic `z` direction.
At steady state, the measured `w(x)` must satisfy

```
(1/Re) d2w/dx2 = b(theta,s),       w(0)=w(1)=0,
b = -betaT |theta-Tmd0-eos_Tmd_slope*s|^q + betaS*s.
```

Two branches distinguish all EOS parameters: a thermal ramp with `s=0`, and a salinity ramp with
fixed temperature. The discrete curvature of `w` reconstructs the body force pointwise.

**Result after correction.** The thermal and salty measured/BVP amplitudes are respectively
`0.99986457` and `0.99991996`; last-output changes are `6.05e-5` and `3.64e-5`. Reconstructed
thermal parameters are `betaT=1.0002443`, `Tmd0=0.4000146`, `q=2.0004528`. Reconstructed salty
polynomial coefficients are `(-0.3600050, 1.1499418, -0.2499509)`, versus nominal
`(-0.36, 1.15, -0.25)`.

The initial runs produced exactly 0.75 amplitude because `Velocity_add_buoyancy_2_RHS` omitted the
periodic wrap face: with `NZM=4`, only three of four `w` faces received buoyancy. An independent
fixed-pressure-gradient Poiseuille test gave unit response and excluded the momentum/viscous
operator. Component-specific periodic buoyancy endpoints were added to match the normal momentum
RHS; the corrected profiles then converged to unit amplitude.

**Gate:** `PARTIES/testcases/1DNonlinearEOS/`; combined analyzer and before/after plots in
`PARTIES/testcases/periodic_fix_validation/`. Full record:
[STAGE_A_EOS_ICE_PENALIZATION_GATE_REPORT.md](STAGE_A_EOS_ICE_PENALIZATION_GATE_REPORT.md).

### A.3d Intermediate validation — ICE penalization / Brinkman sweep  ✅ PASSED 2026-07-05

**Purpose.** Certify A.1.5 using steady Brinkman-Poiseuille flow through a liquid gap adjacent to
a diffuse ice slab. With melting disabled and equal phase viscosities, the measured profile must
satisfy the discrete BVP

```
(1/Re) d2w/dx2 - (1-F)/darcy_tau w = 2s,
```

so the Darcy term is isolated from the high-viscosity rigidity backup. Sweep the nominal
Brinkman thickness `delta=sqrt((1/Re) darcy_tau)` over 2, 4, and 8 grid cells.

**Corrected full-sweep result:**

| `delta/dx` | measured/BVP amplitude | fitted discrete `delta` error | far-ice leakage |
|---:|---:|---:|---:|
| 2 | 0.99999698 | 2.91e-5 | 4.19e-7 |
| 4 | 0.99997370 | 4.06e-6 | 2.52e-5 |
| 8 | 0.99994732 | 1.89e-6 | 2.59e-4 |

All last-output changes are below `5.65e-5`, and all scaled profile residuals are below
`2.63e-6`. The measured slip exponent is `0.597270`; the matching diffuse discrete BVP gives
`0.597286`. The ideal sharp-interface exponent is not the correct finite-band reference for the
2-dx case. The same periodic buoyancy endpoint correction identified by A.3c removes the original
0.75 absolute-amplitude factor; the Darcy coefficient and decay scaling themselves were already
correct. A corrected 1-vs-4-rank continuation agrees to `1.92e-12` in the worst stored field
(`w`: `1.29e-13`).

**Gate:** `PARTIES/testcases/1DIcePenalization/analyze_corrected_sweep.py`; results in
`corrected_penalization_results.csv`, `corrected_penalization_profiles.csv`, and the
`corrected_penalization_*.png` figures. Full record:
[STAGE_A_EOS_ICE_PENALIZATION_GATE_REPORT.md](STAGE_A_EOS_ICE_PENALIZATION_GATE_REPORT.md).

### A.4 Full 2-D Cartesian Yang production run (the single Stage-A benchmark simulation)

**Entry condition:** ✅ SATISFIED 2026-07-05 — A.3–A.3d all passed and the A.2 full-flag
integrated shakeout (`ECCO_TESTS/Yang_checkout/`) passed. The first production run subsequently
failed the joint A.4 acceptance on melt rate, so Stage A remains open pending the sensitivity and
closure runs below.

> **First production run completed 2026-07-12: overall benchmark FAIL.** The 2-D
> `TWOD_CARTESIAN` 1440² run reached `t=200` in 400,019 iterations (2 d 20 h 32 min on
> 144 ranks) without a solver failure. Salt drift was at most 1.33·10⁻¹⁵, F stayed exactly in
> [0,1], and the maximum projection divergence was 7.59·10⁻¹⁰. **Morphology is a qualified
> PASS:** representative scallop spacings are 0.249–0.296H against the H&T prediction 0.287H,
> with the expected thick upper ice and a late merged upper interval of 0.416H. **The required
> melt curve FAILS:** `t½=179.848` vs 118 and `V(200)/V₀=0.45616` vs 0.32 [**targets void — those are the `ΔSv=0` curve; corrected target `t½≈216`, so this run is 1.20× too FAST, see YANG_REFERENCE_CORRECTION_REPORT.md**]. The error is not a
> constant developed-regime slope deficit: the run is 39 % slow over t=40–100, whereas its
> t=120–200 slope is close to Yang's, so the permanent offset forms mainly over t=40–120.
> The leading test is therefore three `darcy_tau={1e-4,2.5e-4,5e-4}` restarts from
> `Resume_10.h5` at t=20 (not `Resume_50` at t=100), through t=80. Full evidence and next steps:
> [YANG_PRODUCTION_ANALYSIS_REPORT.md](YANG_PRODUCTION_ANALYSIS_REPORT.md).
>
> **Branch outcome 2026-07-13: penalization REFUTED as the cause.** τ = 2.5·10⁻⁴ closes 0.8 %
> of the t = 80 volume gap; the τ = 5·10⁻⁴ (η = dt proxy) partial window overlays the control;
> the restart control reproduces production to 10⁻⁸; solid creep grows ∝ τ with no melt
> response. Per the report's Step-2 decision rule: `darcy_tau = 1e-4` stays, and the
> investigation moves to the Allen–Cahn vs Cahn–Hilliard interface-coupling difference under
> convection (phase-field thickness / mobility / convective source-mapping audit + a short
> matched convective-melt case). Details in the report's "Branch outcome" subsection.

> **Final freshwater discriminator 2026-07-21: `Sm = 0` baseline PASS; salty/fresh ratio
> FAIL.** Anvil job 19359337 (`1440²`, 256 ranks) completed `t = 200` cleanly in 39 h 05 min.
> The true `S ≡ 0` half-melt time is `t½,0 = 60.2866`, about 2.2 % from the Yang value
> `≈59` inferred from figure 4(b) [**inference void — figure 4(a) plots the `Sm=0` curve directly at `t½ = 107.42`, so PARTIES is 1.78× too FAST; the 59 was `0.5 × 118` from the same misidentified curve**]. The paired PARTIES result is nevertheless
> `f̄5/f̄0 = 60.2866/179.8484 = 0.3352`, versus Yang `≈0.5`. This directly validates the
> freshwater Yang geometry and localizes the remaining failure to salt-coupled film/layer
> transport or reference uncertainty. The next action is a matched-interface transport audit
> plus the author data request, followed only if needed by short salt-side sensitivities; do
> not repeat Darcy/CH-mobility branches or launch another full run yet. Full record:
> [YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md).

> **Historical allocation-free transport audit status 2026-07-21 (superseded below):** the reusable matched-volume analyzer
> is complete and the Sm=0 states at requested `V/V0={0.8,0.7,0.6,0.5}` have been processed.
> Reconstructed instantaneous melt and finite-difference volume loss agree within 0.12--0.60 %,
> establishing a trustworthy freshwater transport baseline. The Sm=5 raw files remain on the
> separate `/bigscratch` system; run the analyzer there and return its small CSVs rather than
> copying approximately 60 GB or starting a simulation. The next possible job is conditional:
> at most one Sm=5 `Sc=500` restart from `t=20` to `t=60`, and only if the paired audit diagnoses
> a marginal/leaky salt film. Require 20--25 % closure of the existing `t=60` volume gap to
> continue and stop below 5 %. No full or multi-parameter run is currently authorized by the
> evidence. Script, command, outputs, and thresholds are recorded in section 7 of the Sm=0 report.

> **Closure 2026-07-21: salt-resolution and salt-only corrections REJECTED.** The
> returned `Sm=5` audit found all salinity profiles resolved (median recovery length
> 27.7--98.4 cells), negligible salt overlap/leakage, and no starved/refreezing rows. Reduced heat
> flux and interface superheat explain the low melt rate, so the exact decision rule rejects the
> `Sc=500` branch. A subsequent four-rank Le=100 gate mapped the AFiD/Yang finite-interface salt
> equation onto PARTIES. The stable implementation conserved well, but **worsened** both the
> analytic liquid-salinity profile and front constant from production-like and smooth similarity
> starts. Therefore do not run its refined gate, and do not enable `yang_salt_transport` in
> production. Full evidence: [YANG_SALT_OPERATOR_GATE_REPORT.md](YANG_SALT_OPERATOR_GATE_REPORT.md).

> **Controlling closure 2026-07-27 (supersedes the 2026-07-21 route decision).** The action item
> attached to the closure above — "switch to AFiD-MuRPhFi, or port Allen–Cahn before any further
> PARTIES run" — is **withdrawn**; see the ABSOLUTE GOAL block at the top of this document. Two
> hypotheses were never tested and are both cheaper than an Allen–Cahn port:
> **(i)** PARTIES has never been shown to be *grid/interface-width converged* in the salty case —
> `1440²` was adopted because it matches Yang's refined mesh, not because a PARTIES self-convergence
> study demanded it; **(ii)** the reference case sits at/near the **minimum** of the paper's
> non-monotonic `f̄(Sm)` curve (`EXPECTED_RESULTS.md`: minimum near `Sm ≈ 6` for `ΔSv = 5`), i.e.
> the single most delicate point in the paper, and PARTIES has only ever been compared at that one
> point plus the `Sm = 0` endpoint (which it **passes** to 2.2 %). The controlling plan is now
> **§A.5**: settle convergence first, then test the *shape* of PARTIES' own `f̄(Sm)` curve against
> Yang figure 4, then — and only then — consider interface-model form.

**Reference case:** `H = 5 cm`, `ΔT = 20 K` (water 20 °C, ice 0 °C), `Sm = 5 g/kg`,
`ΔSv = 5 g/kg`, `RaT = 10⁷`, `Pr = 10`, `Sc = 1000` (Le = 100), `St_Yang = 4`.

**Grid — why 1440 in x,y.** Yang's Figure 1 grid study is explicitly for this `Sm=5`,
`ΔSv=5` case: it chooses AFiD base resolution 288 for velocity/temperature while fixing the
salinity/phase grid at five-times refinement, `5 × 288 = 1440`. The separately stated
432²/2880² grid is for `RaT=10⁸`, `RaS=2×10¹⁰`, not this reference case. PARTIES is single-grid,
so it uses the full refined `1440²` resolution for every field. The `TWOD_CARTESIAN` build stores
one periodic bookkeeping slab, `1440 × 1440 × 1 ≈ 2.1 M` cells.

**Physical geometry:** this is the paper's 2-D benchmark represented in PARTIES by a one-cell
periodic storage slab. It should be compared with the paper's 2-D data, not with the wall-bounded
3-D case (`Γy=0.5`).

**Numerical-method differences that must remain visible in the comparison:** Yang et al. set
their phase-field thickness parameter to the grid spacing, use `η=dt` for penalization, and add
direct zero-velocity forcing for solid fraction above 0.9. PARTIES uses its independently
validated CH mapping (`Cn=0.75Δx`, `Pe_CH=0.9/Cn`) and fixed implicit
`darcy_tau=1e-4`, without that direct-forcing cutoff. The physical control parameters match;
the phase-field and rigidity implementations are not asserted to be algebraically identical.

**`Boundary.h`:** `TWOD_CARTESIAN`; `x` = no-slip walls (left wall + wall behind the ice),
`y` = no-slip walls (top/bottom), `z` = **periodic** storage slab. `#undef SURFACE_TENSION`. Define `VOF_DIFFUSE`,
`CONC`, `BOUSSINESQ`, `EOS_NONLINEAR`, `CONC_VOF_PHASEWEIGHTED`, `PHASE_CHANGE`,
`ICE_PENALIZATION`, `USE_HYPRE`.

> **Inputs re-audited against the published JFM paper 2026-07-05**
> ([doi:10.1017/jfm.2023.582](https://doi.org/10.1017/jfm.2023.582), plus
> `docs/md_files/Yang2023.pdf`): geometry (Γ = 1, ice slab 0.1H, Ti = 0 °C, all walls
> no-slip/no-flux), EOS eq. (2.8) coefficients (Cb = 0.011, b0 = 0.77, T0 = 4 °C, cS = −0.25 ⇒
> βS = 2b0·Sm/(Cb·ΔT²) = 1.75, Tmd0 = 0.2, Tmd_slope = −0.0625), Pr = 10, Sc = 1000,
> St = L/(cpΔT) = 4 ⇒ `stefan = 0.25`, Re = √(RaT/Pr) = 1000, stratification s ∈ [0.5, 1.5].
> The BC arrays below were corrected to the code's A·∂c/∂n + B·c = C convention (the earlier
> draft had A/B swapped — Dirichlet s = 0 at every wall would have drained the salt field).
> **The ready-to-run as-built input is `PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/`**
> (`parties.inp`, `jobscript.sh`, `EXPECTED_RESULTS.md`).

```ini
[geometry]
# Lx=H=1, Ly=H=1, Lz=one-cell TWOD_CARTESIAN storage slab
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0
zmin = 0.0
zmax = 0.0006944444444444445
twod_slab_thickness = 0.0006944444444444445  # 1/1440; explicit bookkeeping thickness

[grid]
NXM = 1440
NYM = 1440
NZM = 1                 # required by TWOD_MODE

[simulation]
# published volume-curve comparison window ends at t/t_ff = 200
time_max = 200.0
# keep the first step below the fine-grid convective CFL; 1e-3 is too large here
default_dt = 1.0e-4
max_dt = 5.0e-4
cfl = 0.3
output_time_interval = 2.0  # ~101 output times; ~60 GB including Data + Resume files

[flow]
Re = 1000.0             # = sqrt(RaT/Pr) = sqrt(1e7/10)
vel_init_type = 0       # quiescent

[vof]
# water = phase1 (F=1), ice = phase2 (F=0); Boussinesq rho2=rho1
rho1 = 1.0
rho2 = 1.0
mu1 = 1.0
mu2 = 1.0               # constant viscosity in Yang's momentum equation
init_type = 27          # NEW: vertical ice slab, F=0 for x>=0.9
Cn = 5.2083333e-4       # PARTIES-calibrated 0.75*dx; not a literal AFiD epsilon mapping
Pe_CH = 1728.0          # = 0.9/Cn, the same rule as all passed gates (307.2/614.4/153.6)
darcy_tau = 1.0e-4      # validated fixed PARTIES value; Yang uses eta=dt + direct forcing

[phase_change]
stefan = 0.25           # = cp*dT/L = 1/St_Yang ; St_Yang=4
T_melt = 0.0            # theta of ice melting at S=0
liquidus_slope = 0.014  # = m_freeze*Sm/dT; published paper gives m=0.056 C/(g/kg)
melt_band_eps = 5.2083333e-4   # = Cn

[eos]
# Roquet quadratic; nondim s=S/Sm, theta=(T-Ti)/dT
eos_q = 2.0
eos_betaT = 1.0         # quadratic thermal coeff in free-fall units
eos_betaS = 1.75        # = LambdaT = RaS/RaT = 2*b0*Sm/(Cb*dT^2)
eos_Tmd0 = 0.2          # (T0-Ti)/dT = 4/20
eos_Tmd_slope = -0.0625 # cS*Sm/dT (density-max shift; NOT the liquidus)

[conc]
# field 0 = T, field 1 = S
NConc = 2
Pe = {1e4, 1e6}         # Re*Pr , Re*Sc
richardson = {0.0, 0.0} # linear buoyancy OFF: EOS_NONLINEAR supplies b(theta,s)
kappa_ice_ratio = {1.0, 0.0}  # CONC_VOF_PHASEWEIGHTED per-field: T diffuses in ice, S masked
conc_init_type = {30, 31}  # NEW: T step ; S linear-in-y, 0 in ice
cbd0 = 1.0              # theta_water = (Tw-Ti)/dT = 1
theta_ice = 0.0         # ice at the fresh melting temperature (Yang: Ti = 0 C)
cbd2 = 0.5              # s_top = (Sm - dSv/2)/Sm = 2.5/5
cbd5 = 1.5              # s_bot = (Sm + dSv/2)/Sm = 7.5/5
# no-flux (Neumann) for BOTH fields on ALL non-periodic walls.
# Code convention (see A.3 correction note): A*d(c)/dn + B*c = C, A multiplies
# the GRADIENT, so no-flux = A=1, B=0, C=0 (as in every passed gate input).
BC_AN = {1, 1}
BC_AS = {1, 1}
BC_AE = {1, 1}
BC_AW = {1, 1}
BC_BN = {0, 0}
BC_BS = {0, 0}
BC_BE = {0, 0}
BC_BW = {0, 0}
BC_CN = {0, 0}
BC_CS = {0, 0}
BC_CE = {0, 0}
BC_CW = {0, 0}
```

**Timestep check:** PARTIES starts from `default_dt`, then recomputes `dt` from `Dtime_cfl` once
`time != 0`. For the Yang production grid, `dx = dy = 1/1440`; with `cfl = 0.3`, the advective
limit is approximately `dt <= 0.3/(1440*(|u|+|v|))`. Once the melt-driven flow reaches O(1)
free-fall velocity, this is O(`1e-4`) in 2-D and smaller if several velocity components are
active. Therefore `default_dt = 1e-4` is the safe production start. Do not use `1e-3` for A.4; it is
a coarse-grid/debug value only. The `max_dt = 5e-4` cap is acceptable only because adaptive CFL will
choose a smaller value when the flow accelerates.

**Crash fix note (2026-07-09):** the first production submission was mistakenly run as a
four-cell thin 3-D slab (`TWOD_CARTESIAN` undefined, `NZM=4`). It failed at the first pressure
projection with `Poisson solver did not converge` because the 3-D pressure solve did not reduce
the divergence below the hard threshold after 10 correction passes. The corrected production
configuration is `TWOD_CARTESIAN` with `NZM=1`; a reduced 128² smoke test passed with
`TWOD_MODE active` and projection divergence O(`1e-15`).

**Paper-grounded numerical targets — CORRECTED 2026-07-27.** The earlier statement here
(*"digitizing the published 2-D curve gives `t_1/2 ≈ 118 t_ff` … and `V(200)/V0 ≈ 0.32`"*) read the
**`ΔSv = 0`** family of figure 4(a) as if it were `ΔSv = 5`. Those two numbers are the
`Sm = 5, ΔSv = 0` case. This run is `Sm = 5, ΔSv = 5`, whose only published datum is
figure 4(b)'s `f̄/f̄₀ = 0.497`, giving **`t½ ≈ 216`** with `f̄₀ = 1/107.42` measured directly off the
`Sm = 0` curve. Corrected table: `EXPECTED_RESULTS.md` §1.1a and
[YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md).

The layering target is unaffected: using the full EOS for the ambient gradient gives
`|dρ/dz| ≈ 72.26 kg m⁻⁴`; (3.1) then predicts `h ≈ 1.43 cm = 0.287H`, or 3–4 layers.

**Pass criterion (the actual benchmark check):** the melt-front layer thickness vs. ambient
density gradient agrees with Huppert & Turner / published Yang Fig. 3(b) (arXiv Fig. 2b), and
the normalized ice volume matches the published value **for the same `ΔSv`** — figure 4(a) curves
only for `ΔSv = 0` cases, figure 4(b) markers for everything else. The 2-D Cartesian run is
compared with the paper's 2-D data (figure 4(a) curves and figure 4(b) *circles*; the squares are
3-D).

**Stage-transition note for IBM:** Stage A is the **non-particle / non-IBM closure study**. It is
considered over only after A.2 passes the code-shakeout checks, A.3 and A.3b pass the 1D Stefan
gates, and
the single A.4 Yang run passes
the melt-front and ice-volume criteria above. Until then, keep `VOF_IBM`, `LAG_PARTICLE_RESOLVED`,
and `PARTICLE_RELEASE` undefined. Once Stage A is closed, freeze/tag that non-particle state; the
IBM situation then starts in Stage B by enabling `VOF_IBM` and rechecking the ice penalization with
the sediment mask `φ_s,f = (1−F_f)(1−C_{S,f})`, so the Darcy term damps ice but does not add
artificial drag inside the resolved rock particle.

**Reusing for other Yang points:** `eos_betaS`, `eos_Tmd_slope`, and `liquidus_slope` all fold in
`Sm`/`ΔT`. For a different `(Sm, ΔT)` recompute `ΛT = 2 b0 Sm/(Cb ΔT²)`, `eos_Tmd_slope = cS·Sm/ΔT`,
`liquidus_slope = m·Sm/ΔT`, and `Re = √(RaT/Pr)`.

---

### A.5 PARTIES-only closure plan for the Yang melt-rate miss (controlling, 2026-07-27)

This section replaces the withdrawn "run AFiD / port Allen–Cahn" route. It is ordered by
**evidence per core-hour**: free work first, then the cheapest simulation that can change the
verdict, and the interface-model rewrite last — because it is the most expensive hypothesis and
the *least* supported by current evidence.

> **REFERENCE CORRECTION 2026-07-27 — read this first.** §A.5.1's free reference audit was
> executed and it **invalidated the campaign's comparison targets**. Yang figure 3(a)
> (= JFM figure 4(a)), the curve this project digitized, is the **`ΔSv = 0`** series, not
> `ΔSv = 5`. The PARTIES production case is `Sm=5, ΔSv=5`, so it was judged against a different
> physical case, and the freshwater reference was *derived* from that same broken chain.
> Corrected, directly digitized targets: `t½(Sm=0) = 107.42` (was: 59 inferred) and
> `t½(Sm=5, ΔSv=5) ≈ 216` (was: 118). **PARTIES melts too fast in both cases — 1.78× in
> freshwater, 1.20× in the salty case — not "52 % too slow".** Full provenance, three independent
> proofs of the series identification, and the corrected verdict table:
> [YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md). The step ordering
> below is updated accordingly: the **freshwater** case is now the primary target, and it is by
> far the cheapest one to attack because `S ≡ 0` removes the `Sc = 1000` mesh requirement
> entirely.

#### A.5.-1 SYNTHESIS 2026-07-28 — compare at matched EFFECTIVE resolution

Yang's `u,T` sit on **288²** in *both* cases (their 5× refinement covers only `S` and `φ`).
PARTIES freshwater is strongly resolution-sensitive; PARTIES salty is **not** (§A.5.3). Comparing
each code at the resolution that actually sets its heat transport:

| compared at | fresh `t½` | salty `t½` | `f̄₅/f̄₀` |
|---|---:|---:|---:|
| PARTIES 288² | 102.36 | 188.39 | **0.543** |
| PARTIES 1440² | 60.29 | 179.85 | 0.335 |
| **Yang** | **107.42** | **216 ± 11** | **0.497** |

| | fresh | salty | ratio |
|---|---:|---:|---:|
| **at `N=288` (matched)** | **−4.7 %** | **−12.8 %** | **+9.3 %** |
| at `N=1440` (mismatched) | −43.9 % | −16.8 % | −32.6 % |

**The `f̄₅/f̄₀ = 0.335 vs 0.497` failure — treated as the campaign's core unexplained result for
weeks — is largely a resolution-mismatch artifact.** At matched resolution it is **+9.3 %**, not
−32.6 %. The mechanism is now obvious in hindsight: the ratio's *numerator* (freshwater) is
resolution-sensitive while its *denominator* (salty) is not, so computing the ratio at 1440² and
comparing it to a ratio computed at 288² is not a like-for-like comparison.

**Standing at matched effective resolution:** morphology **PASS** (0.249–0.296H vs H&T 0.287H);
freshwater `t½` **−4.7 %** (vs ≈2.8 % digitization uncertainty); ratio **+9.3 %**; salty `t½`
**−12.8 %**, i.e. **2.6σ** against a target that is itself *derived* (`107.42 / 0.497`) and carries
**±11 (4.9 %)**. The residual is one number, at ~2.6σ, against the weakest of the three references.

**Caveat that keeps this honest:** PARTIES at 288² carries a diffuse band **5.3× wider** than
Yang's (§G.8), so both codes under-deliver interfacial heat at that effective resolution but
plausibly for different reasons. The agreement may be numerically right and mechanistically
coincidental. Step C was designed to separate that and failed on the enthalpy leak (§A.5.5 item 3);
until it is redone, treat this as strong circumstantial evidence, not proof.

#### A.5.0 What is already excluded (do not re-test)

| Suspect | Status | Evidence |
|---|---|---|
| Phase-change kernel (G.2–G.3) | **EXCLUDED** | A.3 1-D Stefan `λ` to 1.7 %; in-situ band audit: band superheat = kinetic requirement, `α_eff = α_input`, early melt on the reservoir similarity law (`λ = 0.122` predicted, 0.127 measured) |
| Salt–liquidus coupling | **EXCLUDED** | A.3b 1-D binary Stefan gate; Le=100 gate `λ` error 0.57 % |
| Darcy / ice rigidity | **EXCLUDED** | `τ = {1e-4, 2.5e-4, 5e-4}` branches close 0.8 % of the `t=80` gap; solid creep ∝ τ with **zero** melt response |
| CH mobility `Pe_CH` | **EXCLUDED** | `Pe_CH = 8640` (5×, = Yang's `D` scale) moves `V` by `+1.4e-4` at `t=40` |
| Generic convective transport | **EXCLUDED** | Flat-wall 2-D RB `Nu = 8.27 ± 0.26` at `Ra=1e6` (at/above literature; `Nu_bot=Nu_top=Nu_vol` to 0.1 %) |
| Salt-free convective **melting** | **EXCLUDED** | Favier melting-RB `St=0.1`: plateau ratio 0.98–1.03 of target, diffusive phase 1.1 % off exact 1-D Stefan, budget `Pe·dE/dt/(q_b−q_t) = 1.0006 ± 0.013` |
| Yang geometry + EOS + penalization, freshwater | ~~EXCLUDED~~ **REOPENED 2026-07-27** | the "+2.2 %" pass rested on a misidentified reference curve. Against the directly digitized `ΔSv=0` curve, `t½ = 60.29` vs **107.42** — **1.78× too fast**, now the campaign's largest single discrepancy |
| Unresolved / leaky salt film | **EXCLUDED** | matched-state audit: `δ_S,90` median 27.7–98.4 cells, none `< 4dx`, leakage index 1.25e-4–1.92e-4, no starved/refreezing rows |
| Salt operator form (Yang finite-interface) | **EXCLUDED — and PARTIES won** | Le=100 gate: mapped operator 1.5–1.99× worse on `λ` and 2–6× worse on liquid-`S` `L∞` from both a production-like and an exact-similarity start |

**Everything tested against an *exact analytic* reference still passes** — the 1-D gates, the
melting-RB budgets, the flat-wall Nusselt. What fails, against the corrected digitized reference,
is the **melt rate in the Yang configuration itself**, in both the freshwater and the salty case:

| | PARTIES | corrected Yang target | |
|---|---:|---:|---|
| `t½ (Sm=0, ΔSv=0)` | 60.287 | **107.42** | 1.78× too fast |
| `t½ (Sm=5, ΔSv=5)` | 179.848 | **≈216** | 1.20× too fast |
| `f̄₅/f̄₀` | 0.3352 | **0.497** | 0.67 of target |
| melt-front layer scale | 0.249–0.296 H | 0.287 H | **PASS** |

The morphology pass is important and unaffected: it independently confirms the ambient density
gradient, hence the stratification magnitude, the `ΔSv` sense, and the EOS are set up correctly.
The `ΔSv` sense is separately confirmed from the paper's initial condition (§2 of the correction
report).

Note what the corrected numbers do to the *shape* of the problem: the salty case is now the
**closer** of the two, and the freshwater case — which contains no salt physics at all — is the
outlier. Any explanation that lives in the salt coupling cannot account for the freshwater miss.

#### A.5.1 Step 0 — free work ✅ DONE 2026-07-27 (items 1–3)

Items 1–3 are complete and are reported in
[YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md). Summary of what the
paper itself settled, at zero allocation and with no author data:

- **`f̄ = 1/t½`** — verbatim from `main.tex`. The project's convention was right.
- **`ΔSv` sense** — the initial condition `S = S_bot + (S_top − S_bot) z/H` plus "stronger stable
  stratification" fixes `S_bot > S_top`. PARTIES' `cbd2 = 0.5`, `cbd5 = 1.5` was right.
- **Figure 3(a) is 2-D** (its ratios match panel (b)'s circles, not its squares), so the
  `TWOD_CARTESIAN` comparison is dimensionally correct.
- **Figure 3(a) is the `ΔSv = 0` series** — the finding that invalidated the targets. Three
  independent proofs in §3 of the report.
- The arXiv source is archived at `docs/papers/yang2023_arxiv_source/`; the digitizer is
  `Yang_production/digitize_yang_fig3.py`; corrected data are `yang_fig3a_reference.csv` (all
  four `Sm` at `ΔSv=0`) and `yang_fig3b_reference.csv`.
- **The Lohse data request is no longer blocking** and need not be sent. All four of its questions
  are answered. Note the letter's own text has an internal inconsistency (it describes figure 3(a)
  as `ΔSv=5`), so author data would still be *nice* — but nothing waits on it.

Remaining Step-0 item — **instrumented 2026-07-27, measurement still owed:**

4. **Audit the conservative clip in 2-D.** `VOF_DIFFUSE_step` calls
   `diffuse_redistribute_liquid_mass_delta`, which is **globally** MPI-collective: clipped melt
   mass is redeposited anywhere on the interface with `F(1−F)` headroom, not where it was removed.
   In 1-D this is provably a no-op; in the salty 2-D run the row-melt CV is 0.35–0.64, so a
   non-zero `clip_delta` would *teleport* melt from fast rows to slow rows and flatten exactly the
   scallop dynamics that set heat delivery.

   **Instrumentation is in place**: `VOF_DIFFUSE_CLIP_AUDIT` in
   `src/Eulerian/VOF_DIFFUSE.c` (top of file, commented out by default; enable with
   `-DVOF_DIFFUSE_CLIP_AUDIT=100` or by uncommenting). It prints per-step `clip_delta` and the
   redistribution `residual`, plus running `Σ|clip|`, `Σ|residual|` and `Σ melt·dt`, and their
   ratio. Log-only: with the macro undefined the block differs from the original by one unused
   local, so the production build is unchanged. Both paths syntax-check.

   **✅ MEASURED 2026-07-28** (`ECCO_TESTS/StageA_Yang/Yang_clip_audit/`, `parties.clipaudit`, on-mapping
   control vs the Step-C off-mapping config, both to `t=3`):

   | `Σ|clip| / Σ melt` | ntime=100 | 500 | 1000 | 1500 | enthalpy drift at `t=3` |
   |---|---:|---:|---:|---:|---:|
   | control (`Cn = 0.75Δx`) | 4.5e-3 | 1.6e-2 | 2.1e-2 | **2.4e-2** | **1.8e-7** |
   | off-mapping (Step C) | 9.9e-4 | 4.7e0 | 8.5e0 | **1.33e1** | **8.1e-3** |

   **Hypothesis CONFIRMED.** Off the tuned mapping the clip moves **13× more mass than the total
   melt**, 555× the control, and the enthalpy drift tracks it (4.5 orders). The redistribution
   *residual* is ~1e-22 in both, so `F` mass is placed successfully — the leak is not failed
   redistribution but the fact that the clip/redistribute cycle relocates `F` **after** the
   temperature field has already paid the latent heat for it.

   > **Correction to §G.3.** The claim there that "with `w(F)` the clip stays a no-op" is
   > **inaccurate**. Even on the tuned mapping the clip handles **2.4 % of cumulative melt** and is
   > still growing at `t=3`. What is true is that the conservative redistribution *absorbs* it —
   > drift stays at 1.8e-7. The clip is small and absorbed, not absent. That distinction matters:
   > it is the margin that disappears off-mapping.

   Acceptance was predeclared as `≤ O(1e-6)`; the control misses that by four orders yet is
   demonstrably healthy, so **the threshold was wrong, not the code**. The operative criterion is
   the enthalpy drift the redistribution leaves behind (`≤1e-6`), with the clip ratio as the
   leading indicator. A band-local restore (as in the trapped-cavity work:
   `VOF_DIFFUSE_BAND_MASS_RESTORE`, `diffuse_band_conserve_contact_line`, **not** present in this
   tree) would raise that margin, and is the right fix if any future configuration needs to leave
   the `Cn = 0.75Δx` mapping.

   *Do not rebuild the repository binary while the A.5.2 ladder is in flight*: those runs copied
   `parties` into their own directories, so they are safe, but the ladder's comparability argument
   rests on one binary throughout.

#### A.5.2 Step 1 — freshwater convergence ladder ✅ RUN 2026-07-27 — RESULT BELOW

> **RESULT: the freshwater discrepancy is a RESOLUTION comparison, not a PARTIES defect.**
> Six rungs, five of them run with one binary and identical physics, only the grid and the four
> grid-tied numbers varying:
>
> | `N` | 216 | 288 | 360 | 512 | 720 | 1440 |
> |---|---:|---:|---:|---:|---:|---:|
> | `t½` | 116.31 | **102.36** | 94.47 | 85.59 | 78.89 | 60.29 |
>
> **Yang's digitized `t½ = 107.42` falls between the `N=216` and `N=288` rungs**, i.e. it equals
> PARTIES at an effective `N = 257`. Yang's stated base grid for **velocity and temperature** in
> this case is **288²** — 10.8 % away — and their 5× refinement to 1440² applies only to salinity
> and the phase field, which buys nothing here because the `Sm=0` case **has no salinity at all**.
> PARTIES run at Yang's own base resolution gives `102.36` against their `107.42`, a **4.7 %**
> difference (figure digitization uncertainty on `t½` is ≈ ±3, i.e. 2.8 %).
>
> The ladder is monotone, smooth, and predictive: it forecast `N720` at 79.5 and got 78.89 (0.8 %),
> and forecast `N288` at 101–102 and got 102.36. Extrapolating `t½` in `1/N` to infinite
> resolution gives **`t½(∞) ≈ 63–64`**.
>
> **Mechanism, established at matched melt state (`V/V₀ ≈ 0.75`):** interface arclength is
> grid-independent to 0.15 % and the front is nearly flat, so this is *not* exposed area — it is
> heat delivered per unit area. The liquid thermal layer `δ_T` is well resolved on every rung
> (18–42 cells). What changes is the diffuse band as a **fraction** of that layer —
> 16.7 % → 11.1 % → 7.1 % at `N` = 360/512/720 — because `Cn = 0.75Δx` keeps the band at 3 cells
> while the physics stays fixed. A finite band smears the interfacial temperature gradient and
> suppresses heat delivery, so coarse grids melt slower and look deceptively closer to Yang.
>
> **The `1440²` point is REAL — provenance PASSED, and the ladder's extrapolation is what fails.**
> An earlier draft of this block claimed the `1440²` point was "inconsistent with its own ladder"
> because it sits below the extrapolated `t½(∞) ≈ 63–64`, "which a monotone convergence sequence
> cannot do", and suspected the different binary (`aac878bf` vs the current `6a8fa6eb`). **That
> reasoning was backwards.** The provenance re-run (`ECCO_TESTS/StageA_Yang/Yang_fresh_provenance/`, 1440²,
> current binary, 128 ranks vs the original 256) reproduces the completed run to
> **1.7e-12 at `t=2`**, 1.1e-5 at `t=4` and 5.8e-5 at `t=6` — against a discriminating size of
> 7.2e-3 at `t=6`, a **125× margin**. The `t=2` agreement is effectively bitwise, proving the code
> paths are identical; the later growth is round-off from a different rank count (different
> summation order) amplifying in a convecting flow. Binary **and** rank-count independence are
> therefore both confirmed, which also clears the lineage under the salty production run
> (2026-07-12), since it predates the current binary too.
>
> **What actually fails is the extrapolation.** Successive slopes `dt½/d(1/N)`:
>
> | rungs | 216→288 | 288→360 | 360→512 | 512→720 | **720→1440** |
> |---|---:|---:|---:|---:|---:|
> | slope | 12058 | 11369 | 10764 | 11879 | **26800** |
>
> The coarse rungs are consistent at ≈11 000–12 000; the last interval jumps by 2.3×. A low-order
> polynomial fitted over 216–720 was never entitled to predict 1440, a factor of two beyond its
> range and in a regime where the sensitivity was visibly growing. **Do not quote `t½(∞) ≈ 63–64`
> as a converged limit** — the ladder is not converged anywhere in the sampled range, including at
> 1440². A plausible physical reading of the slope jump is that the sidewall boundary layer goes
> unsteady between `N=720` and `N=1440` and the transition to time-dependent convection sharply
> raises heat transport; that is a hypothesis, not a finding.
>
> Consequence for the headline: **strengthened, not weakened.** PARTIES' converged answer is
> `≤ 60.29`, i.e. *further* from Yang, and the 288²-bracket comparison uses only rungs inside the
> well-behaved range. What we may no longer claim is a specific converged `t½`.
>
> **⚠ Mechanism vs number.** PARTIES at `N≈258` has a fat band (≈23 % of `δ_T`) whereas Yang keeps
> a sharp phase field (`ε = 1/1440`) and coarsens only `u` and `T`. Both under-deliver interfacial
> heat at that effective resolution, but plausibly for different reasons. The `t½` agreement could
> therefore be numerically right and mechanistically coincidental — which is exactly what §A.5.5
> variant 3 (freeze `Cn`, vary `Δx`) now has to separate. Treat the conclusion as strong
> circumstantial evidence, not proof, until that runs.
>
> Data: `ECCO_TESTS/StageA_Yang/Yang_fresh_convergence/` (`make_inputs.py`, `analyze_ladder.py`,
> `ladder_summary.csv`, `README.md`).

**Original plan text follows.**

**Why this is now first.** The freshwater miss (1.78×) is larger than the salty one (1.20×) and
contains **no salt physics**, so it is both the bigger signal and the simpler system. Critically,
`S ≡ 0` removes the `Sc = 1000` mesh requirement altogether — the `1440²` mesh was demanded by the
salt field, and the freshwater case only has to resolve `Pe_T = 1e4`. This turns the most
expensive question in the campaign into the cheapest one. PARTIES has **never** run a convergence
study of this case; the single completed `Sm=0` run inherited `1440²` from the production input.

- **Decks:** `ECCO_TESTS/StageA_Yang/Yang_fresh_convergence/{N360,N512,N720}/` — generated by
  `make_inputs.py`, each with `parties.inp` + `jobscript.sh`. Every physical input (`Re`, `Pe`,
  EOS, `stefan`, `liquidus_slope`, `darcy_tau`, BCs, geometry, `vof_slab_x0`) is byte-identical to
  the completed `Sm=0` run; **only** the grid and the four grid-tied numbers change:
  `Cn = melt_band_eps = 0.75Δx`, `Pe_CH = 0.9/Cn`, `zmax = twod_slab_thickness = Δx` (full double
  precision — `GRID_UNIFORM` aborts otherwise), and `dt ∝ Δx`. `NConc = 2` with `s ≡ 0` is
  retained rather than `NConc = 1` so the ladder is exactly comparable to the existing `1440²`
  point.
- **Because `Cn = 0.75Δx`, refining the mesh also sharpens the diffuse band** — this is a
  diffuse-interface convergence study, not merely a mesh study.

| `N` | `Cn` | `Pe_CH` | `dt₀` | ranks | `t→` | est. cost |
|---:|---:|---:|---:|---:|---:|---:|
| 360 | 2.0833333e-3 | 432.0 | 4.0e-4 | 36 | 140 | ≈ 110 core-h |
| 512 | 1.46484375e-3 | 614.4 | 2.8e-4 | 64 | 140 | ≈ 315 core-h |
| 720 | 1.0416667e-3 | 864.0 | 2.0e-4 | 144 | 140 | ≈ 875 core-h |
| 1440 | 5.2083333e-4 | 1728.0 | 1.0e-4 | — | 200 | **already run** (`t½ = 60.287`) |

≈ 1 300 core-h for the whole ladder — about one eighth of a single salty production run.
`t → 140` brackets both `t½` values (PARTIES 60.3, Yang 107.4).

- **Precommitted decision rule** on `t½(N)`:
  - **`t½` flat across the ladder** (spread < 5 %, i.e. converged at 360² already) ⇒ resolution is
    not the cause; the 1.78× is a **model or reference** discrepancy. Go to A.5.5 model variants,
    and treat §6 of the correction report as live — i.e. do not assume PARTIES is the party in
    error. The cheap ladder then becomes the sandbox for every variant.
  - **`t½` rising with `N` toward ≈107** ⇒ the completed `1440²` run is under-resolved and the
    remedy is refinement. Richardson-extrapolate and, only then, decide whether a `2880²`
    freshwater confirmation is worth ≈ 7 000 core-h.
  - **`t½` falling with `N`** ⇒ refinement moves *away* from Yang; that is a strong result and
    goes straight to A.5.5.
- **Also record** at matched `V/V₀` on every grid: interface superheat, liquid-side heat flux,
  `max_speed`, and the mean interface length (scallop area). The last one matters — a scalloped
  front has more area than a flat one, and interface area is the most likely route by which
  resolution changes the melt rate.

#### A.5.3 Step 2 — salty resolution test ✅ RUN 2026-07-28 — THE FRESHWATER STORY DOES **NOT** TRANSFER

> **RESULT: the salty melt rate is resolution-CONVERGED, so the remaining 1.20× gap is real and is
> NOT a resolution artifact.** `Sm=5, ΔSv=5` at `N` = 288 and 432 (`Yang_salty_resolution/`),
> against the completed 1440²:
>
> | `N` | 288 | 432 | 1440 | spread | Yang |
> |---|---:|---:|---:|---:|---:|
> | **salty `t½`** | 188.39 | 178.28 | 179.85 | **5.4 %** | ~216 |
> | freshwater `t½` | 102.36 | — | 60.29 | **70 %** | 107.42 |
>
> Over the *same* 5× grid range the freshwater melt rate changes by 70 % and the salty one by 5 %
> (and non-monotonically, i.e. scatter rather than trend). The naive prediction from freshwater
> scaling — `179.85 × 1.698 = 305` — is decisively refuted.
>
> **Health verified at `t½`, not max-over-run** (the distinction that mattered for the ladder, and
> the one my first analyser pass got wrong): enthalpy drift **7.6e-6** (`N288`) and **8.1e-6**
> (`N432`) at `t½`; `∫s` conserved to **4.9e-16**; negative-salinity cells **0.000 %** of the
> liquid at `t½` (0.004 % max, late). The `s_min = −0.055` / `s_max = 1.543` extremes are paired
> local dispersion from central differencing at `Pe_S = 1e6` on a marginal mesh — real, but
> confined to a handful of cells and appearing long after `t½`. The `t½` values stand.
>
> **Why the two cases differ (physical reading).** `ΔSv = 5` stratification suppresses convection
> strongly — that is the paper's own mechanism for the melt-rate minimum. The salty flow is
> therefore laminar enough to be resolved already at 288², whereas the freshwater case has vigorous
> convection whose heat delivery is resolution-limited. Two regimes, two answers.
>
> **Consequence — this SPLITS the campaign cleanly:**
> - **Freshwater: CLOSED.** The 1.78× was a comparison of PARTIES at 1440² against Yang at an
>   effective 288² (§A.5.2).
> - **Salty: OPEN, and now sharply constrained.** PARTIES gives ≈180 at *any* resolution ≥288²;
>   Yang gives ≈216. Resolution is excluded, so this is a model-form or reference difference.
>
> **The new constraint is strong and should drive §A.5.5:** whatever explains the salty gap must
> act **only when salt is present**, because the freshwater case now agrees at matched resolution.
> That promotes **§A.5.5 variant 1 (liquid-referenced salinity `s/(F+δ)` in the liquidus and EOS)**
> to the top candidate — it acts only where `0 < F < 1` *and* salt exists, i.e. exactly and only in
> the salty film. Candidates that would also perturb the freshwater case are now excluded by
> construction.

**Original (pre-result) plan text follows.**

#### A.5.3b Step 2 (original design) — salty grid / interface-width convergence pair

**Demoted below the freshwater ladder** by the reference correction, but still the right test on
the salt side, and A.5.2 will already have calibrated how PARTIES' melt rate responds to
resolution in this geometry. Run it only if A.5.2 shows a resolution response, or if the salty
`f̄₅/f̄₀ = 0.335` vs `0.497` ratio survives whatever A.5.2/A.5.5 establish.

`1440²` was chosen to match Yang's refined mesh, not from a PARTIES convergence study. Because
`Cn = 0.75Δx`, refining the grid also sharpens the diffuse band, so a grid pair *is* the
diffuse-interface convergence test — and on the salt side there is a mechanism that hurts the
salty case much more than the fresh one: the compositional film is `Le^{1/3} ≈ 4.6×` thinner than
the thermal layer, so a 4–6 cell interface band is a few percent of the fresh thermal layer but a
large fraction of the salty film.

- **Run:** `Sm=5` reference input, unchanged except `NXM=NYM=720`, `zmax`/`twod_slab_thickness`
  `= 1/720` (**full double precision — `GRID_UNIFORM` aborts otherwise**), `Cn = melt_band_eps =
  0.75/720 = 1.0416667e-3`, `Pe_CH = 0.9/Cn = 864`. From `t=0` to `t=60`. Everything else fixed.
- **Cost:** ≈ 375 core-h (4× fewer cells, ~2× larger CFL-limited `dt` ⇒ ≈ 1/8 of the ≈ 3 000 core-h
  that `1440²` needs to reach `t=60`). The `1440²` half of the pair **already exists** — free.
- **Comparison point:** `V/V₀` at `t=60`. PARTIES `1440²` = `0.77427`; Yang = `0.69664`; gap
  `= 0.07763`.
- **Precommitted decision rule** — let `Δ = V₇₂₀(60) − V₁₄₄₀(60)`:
  - `|Δ| < 0.008` (< 10 % of the gap) ⇒ **converged.** Interface width is not the cause; the CH
    formulation is exonerated a third time. Adopt `720²` as the validated workhorse grid and go
    straight to A.5.3 at 1/8 cost.
  - `Δ > 0.008` (coarse melts *slower*, i.e. refinement melts faster) ⇒ **live mechanism.**
    Richardson-extrapolate `V∞ = V₁₄₄₀ − Δ/(2^p − 1)` reporting both `p=1` and `p=2`. If
    `V∞ ≤ 0.7549` (≥ 25 % of the gap closed) authorize **one** conditional `2880²` short window
    `t = 20 → 40` seeded by interpolation from `Resume_10.h5` (≈ 8 000 core-h — the only expensive
    run in this plan, and only with this trigger).
  - `Δ < −0.008` (refinement melts *slower*) ⇒ grid is not the cause and refining makes it worse;
    skip to A.5.3 and A.5.4.
- **Also record** at `t=60`, using `analyze_yang_matched_transport.py` on both grids: interface
  superheat, liquid-side heat-flux proxy, `δ_S,90/dx`, and row-melt CV. Those tell you *how*
  resolution acts even if `V(60)` barely moves.

#### A.5.4 Step 3 — reproduce the *shape* of PARTIES' own `f̄(Sm)` curve

The reference correction hands this step a much better target than it had. Figure 3(a) is the
`ΔSv = 0` series, so PARTIES now has **four fully digitized `V(t)/V₀` curves** to compare against
curve-for-curve — not a single ratio point. And PARTIES has run **none** of them except `Sm=0`:
the production case is `ΔSv = 5`, which appears in the paper only as one marker in panel (b).

`ΔSv = 0` means uniform initial salinity, so `cbd2 = cbd5 = 1.0` (`s = S/Sm ≡ 1`). Everything
else re-scales with `Sm` exactly as before (`βS = 2b₀Sm/(CbΔT²)`, `eos_Tmd_slope = cS·Sm/ΔT`,
`liquidus_slope = m·Sm/ΔT`); `eos_q = 2`, `eos_betaT = 1`, `eos_Tmd0 = 0.2`, `stefan = 0.25`,
`Re = 1000`, `Pe = {1e4, 1e6}` are unchanged.

| `Sm` | `eos_betaS` | `eos_Tmd_slope` | `liquidus_slope` | `cbd2`,`cbd5` (`ΔSv=0`) | Yang `t½` | Yang `V(200)/V₀` |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | — | — | — | 0, 0 | **107.42** | 0.2483 |
| 5 | 1.75 | −0.0625 | 0.014 | 1.0, 1.0 | **118.42** | 0.3142 |
| 10 | 3.50 | −0.1250 | 0.028 | 1.0, 1.0 | **112.21** | 0.2837 |
| 15 | 5.25 | −0.1875 | 0.042 | 1.0, 1.0 | **101.84** | 0.2245 |

For reference, the `ΔSv = 5` production case (`cbd2 = 0.5`, `cbd5 = 1.5`) has target
`f̄/f̄₀ = 0.497` ⇒ `t½ ≈ 216`; it has no published `V(t)` curve.

- **Run first:** `Sm = 5, ΔSv = 0`. One run, a **complete curve** target, and it isolates "add
  salt" from "add stratification" — the production case changes both at once. Grid from A.5.2/A.5.3
  (these cases carry salt, so `Sc = 1000` does constrain the mesh; they are not as cheap as the
  freshwater ladder).
- **Then** `Sm = 10` and `Sm = 15` at `ΔSv = 0` for the curve shape.
- **Reading the result:**
  - PARTIES tracks all four `ΔSv=0` curves but misses `ΔSv=5` ⇒ the defect is in the
    **stratification** response, not the salt response. That is a sharp, publishable localization.
  - PARTIES is uniformly fast at every `Sm` including `Sm=0` by a common factor ⇒ a systematic
    transport or scaling offset in this geometry; A.5.2 will already have said whether it is
    resolution.
  - PARTIES reproduces the **non-monotonic ordering** (`Sm=5` slowest, `Sm=15` fastest, minimum
    near `Sm ≈ 3.5`) even if absolute rates differ ⇒ the physical mechanism of the paper is
    reproduced, which is the scientifically substantive claim, and the absolute offset can be
    reported as such.

#### A.5.5 Step 4 — interface-model variants, tested in the 1-D sandbox first

**Never spend 2-D allocation on a model variant that has not first been run in
`ECCO_TESTS/StageA_Yang/Yang_salt_Le100_gate/`** — it has an exact two-phase similarity solution
(`λ = 0.21713826597063368`, `S_int = 1.1640201641615062e-3`), runs on four ranks in ~2 min
(≈ 0.15 core-h), and it already rejected one variant cleanly. Candidate variants, in order:

1. **Liquid-referenced salinity in the band — ✅ IMPLEMENTED AND PASSED THE 1-D GATE 2026-07-28.**
   The liquidus and the EOS used the volume-averaged `s`, but the physical driver is the *liquid*
   salinity `s/(F+δ)`. It acts **only** where `0 < F < 1` **and** salt is present — exactly the
   selectivity the §A.5.3 freshwater/salty split now demands.

   Implemented as runtime flag `[phase_change] liquid_referenced_salinity` (default **0**), in
   `VOF_DIFFUSE_compute_melt_rate` (liquidus) and the `EOS_NONLINEAR` block of
   `Velocity_add_buoyancy_2_RHS` (buoyancy).

   **Le=100 gate result** (`Yang_salt_Le100_gate/`, vs the exact two-phase similarity solution):

   | branch | fitted `λ` | rel. error | liquid-`S` `L∞` | salt drift | max speed |
   |---|---:|---:|---:|---:|---:|
   | flag off (regression) | 0.21844955 | 6.039e-3 | 1.441e-2 | 2.764e-3 | 0 |
   | **flag on (variant 1)** | 0.21749661 | **1.650e-3** | 1.457e-2 | 2.764e-3 | 0 |
   | *mapped Yang salt operator (rejected 07-21)* | *0.21908046* | *8.944e-3* | *8.773e-2* | — | — |

   **The front constant is 3.7× more accurate**; the salinity profile is 1.1 % worse; salt drift
   and quiescence unchanged. Gate verdict: **RUN REFINED GATE**. Note the contrast with the
   rejected Yang salt-operator graft, which made `λ` *worse* — this variant goes the other way.

   **Regression verified:** the rebuilt binary on the original input reproduces the recorded
   2026-07-21 legacy numbers to every printed digit (`λ = 0.21844955`, `L∞ = 1.441e-2`,
   drift `2.764e-3`), so flag-off is behaviourally unchanged.

   **⚠ What this does NOT show.** The gate is quiescent with `eos_betaT = eos_betaS = 0`, so it
   tests **only the liquidus half**. That half is small in production (`Λ* = 0.014`) — the gate
   amplifies it 36× by using `liquidus_slope = 0.5`. The **EOS half is the dynamically important
   one** (`βS = 1.75`) and is completely untested here. Next step is 2-D at `432²`, which §A.5.3
   showed is resolution-converged for the salty case, so ≈470 core-h instead of ≈10 000.

   **Build note (cost a failed job):** do **not** pass `cflags=` on the `make` command line — a
   command-line assignment *overrides* the makefile's target-specific `cflags +=` and silently
   drops the optimisation flags, producing a binary that segfaults. Edit the source or use a
   separate build; `make clean` after any such accident.

   **2-D 432² EOS-half confirmation — ✅ RUN 2026-08-26 (`YangSalty_07282026/N432_sliq1/`) —
   RESULT: NULL, EOS half excluded.**

   > `t_half = 177.99` (flag on) vs `178.28` (flag off, `Yang_salty_resolution/N432`) — a 0.16 %
   > difference, and the wrong direction to close the gap to Yang `≈216.76`. Health checked at
   > `t_half`: enthalpy drift `6.5e-6` (matches the flag-off health, `7.6–8.1e-6`), salinity bounds
   > `s ∈ [-0.0000, 1.49]`, sane. Job ran `t=75.15→280` on the existing checkpoint (originally
   > started 2026-07-28, paused, resumed 2026-08-26); COMPLETED, exit 0.
   >
   > **The EOS/buoyancy half has no measurable effect on the 2-D melt rate.** Combined with the
   > 1-D gate result above (liquidus half: small but real, 3.7× more accurate `λ`, on a term that
   > is small in production), **variant 1 is now fully tested in both halves and fully excluded.**
   > The salty `t½` residual (**−12.8 %, 2.6σ**) stays unexplained; no further leads from §A.5.5
   > remain open.
2. **Band-local mass restore** (from A.5.1.4) if the clip audit demands it.
3. **Interface-width mapping — FIRST ATTEMPT FAILED 2026-07-27, and the failure is informative.**
   Design: mesh fixed at `N=360`, `Cn = melt_band_eps` raised to the `N=216` value
   `3.4722e-3` and `Pe_CH = 0.9/Cn = 259.2`, widening the band from 3 to 5 cells — *better*
   resolved for the CH operator. Expected `t½ ≈ 94.5` (mesh controls) or `≈116.3` (band controls).
   **Got neither: the ice melted out entirely by `t ≈ 10`** (70 % gone by `t=2`, versus 3.7 % in
   the standard `N360`).

   It is **not** a physics result — it is a **conservation failure**. Enthalpy `∫(θ + F/St)` drifts
   `5.2e-3` by `t=2` and `1.1e-2` by `t=8`, against `1.1e-6` for the `N216` run that has the
   *identical* `(Cn, Pe_CH, melt_band_eps)` triple on its own mesh. Four orders of magnitude worse.
   The deck was verified to differ from `N360` in only those three keys plus `time_max`, and the
   melt-law scalars (`m_max = 0.367`, explicit `dt` limit `0.482` ≫ `dt = 4e-4`) are identical to
   `N216`'s.

   **Interpretation:** `Cn = 0.75Δx` with `Pe_CH = 0.9/Cn` is a *tuned combination*, not three free
   parameters. Off that mapping, the melt source outruns the CH relaxation, `F` is driven past 1
   behind the front, and the admissibility clip destroys melt mass whose latent heat the
   temperature field has already paid — precisely the leak documented in §G.3 and the reason
   `w(F) = F(1−F)/(√2 Cn)` was adopted in the first place. `F_max` pinned at exactly 1.0 is
   consistent with that.

   **Consequences.** (a) The "is the knob `Cn` or `Δx`?" question **cannot** be answered by varying
   `Cn` at fixed mesh in this code, so the A.5.2 mechanism-vs-number caveat stays open. (b) The
   `VOF_DIFFUSE_CLIP_AUDIT` instrumentation (§A.5.1 item 4) now has an ideal test case: rebuild with
   it on and re-run this configuration for ~200 steps; it should show a large `clip_delta`, which
   would confirm the mechanism and simultaneously validate the diagnostic. (c) **The A.5.2 ladder
   is unaffected** — every rung sits at ≈1e-6 enthalpy drift at its own `t½`; only this off-mapping
   configuration leaks.

   A corrected design must vary **one** thing at a time. `melt_band_eps` and `Cn` enter separately
   (`coeff ∝ 1/ε`, and `w(F)` is normalised by `Cn`), so the clean test of band width *at fixed
   front speed* is to change `Cn`+`Pe_CH` while holding `melt_band_eps` at the mesh-matched value.
   Do not attempt it before the clip audit explains the present failure.
4. **Allen–Cahn/Hester phase equation — LAST.** Only if A.5.2–A.5.4 and variants 1–3 all come
   back null. It must then be implemented *in PARTIES* (Allen–Cahn
   phase equation, full per-stage phase increment into latent heat and salt in the AFiD substep
   order, liquidus forcing, interface salt-flux term, `δ`, and solid-velocity treatment), verified
   on manufactured arrays, checked for 1-rank vs 4-rank identity and restart history, and passed
   through the Le=100 gate — where it must **beat** the current CH result, which the mapped salt
   operator did not. Budget this as weeks of development, and do not start it before A.5.1–A.5.3
   report.

#### A.5.6 Final production run and Stage-A closure

**The old `t=60` gap thresholds are void** — they were computed against the misidentified
`ΔSv = 0` curve. New gate: a further `1440²`, `t → 200` production run (≈ 10 000 core-h) is
authorized **only** after a step above closes ≥ 25 % of a *corrected* gap, i.e. moves
`t½(Sm=0)` from 60.29 at least 25 % of the way to 107.42 (`t½ ≥ 72.1`), or the equivalent for
whichever case is being targeted, with bounded `F` and `s`, salt drift `< 1e-8`, enthalpy closure,
and coherent layering.

Stage A closes when `t½` (within 10 %), `V(200)/V₀` (within ≈ 0.05), the `0.287H` layer scale, and
all conservation/solver gates pass together **for the case actually being compared** — and the
case must be stated with its `ΔSv`, which the campaign previously did not do.

### ✅ STAGE A CLOSED 2026-07-30

A.5.2–A.5.5 returned null, and the outcome above is the one adopted: **report the discrepancy and
proceed.** Decision taken by the project lead on 2026-07-30.

**Basis.** Every avenue that could close the residual was tested and excluded — mesh resolution
(the salty case is converged), the salt transport operator (ours beats the reference's against an
exact similarity solution), band-localized model corrections (Brinkman drag overlaps the band and
damps them, which constrains *both* codes), Darcy penalization, CH mobility, and salt-film
resolution. What remained — an Allen–Cahn port — costs weeks against a low expected return, and the
one exact reference available indicates our conservative Cahn–Hilliard is the *more* accurate of
the two formulations.

**Amended criterion.** "Within 10 %" is retained for **directly-referenced** quantities and is met:
freshwater `t½` −4.7 %, `V(200)/V₀` −1.4 %, melt-rate ratio +9.3 %, layer scale brackets `0.287H`.
For quantities available only as **derived** references it becomes "documented and bounded": the
salt-stratified `t½` is **−12.8 %, 2.6σ against a target carrying ±4.9 %** (reconstructed as
`107.42 / 0.497` from one curve and one marker). This residual is **unexplained** and is carried
forward in the limitations of the validation summary and the manuscript — closure records it, it
does not dismiss it.

**What this unblocks:** Stage-B *implementation* and the reduced 3-D shakeout. It does **not**
unblock the ECCO production campaign, which stays gated behind that shakeout.

---

## STAGE B — ECCO SEDIMENT-FROM-ICE STUDY

### B.0 Non-invasiveness check (gate on everything below)  — ✅ **BOTH PARTS PASSED 2026-07-30**

Stage B turns on `VOF_IBM` + `LAG_PARTICLE_RESOLVED`. Before building anything on that, the
question that must be answered is: **does enabling those flags perturb the Stage-A physics when no
particles are present?** If it does, every Stage-A validation is invalidated for Stage-B use and
B.1–B.5 are built on sand.

Binary `parties.stageB` (`85e0d33f`, 1 082 184 B vs Stage-A `bc308445` 1 043 336 B). The Stage-A
config and binary were restored immediately after building; `Boundary.h.stageA` is the backup.

The check has **two deliberate parts**, because the Le=100 gate is quiescent (max speed exactly 0)
and so exercises the melt kernel but *not* the flow-coupled IBM paths — the wetting/contact-line
overwrite and the `C_L + C_S <= 1` admissibility trim only bite when there is motion.

| part | test | result |
|---|---|---|
| 1 | Le=100 salt gate vs the recorded Stage-A reference | ✅ **PASS** — `lambda = 0.21844955`, rel. err `6.039e-3`, profile `L_inf 1.441e-2`, salt drift `2.764e-3`: **identical to every printed digit** to the Stage-A value `0.21844955456451656` |
| 2 | freshwater 288² to `t = 10` vs the validated Stage-A run (job 19580986) | ✅ **PASS** — worst relative ice-volume difference **9.2e-9**; see below |

Part 1 shows the two guards that activate under `VOF_IBM` are true no-ops at `C_S = 0`: the melt
source's sediment exclusion (`VOF_DIFFUSE.c:1344`) never fires, and the `C_L + C_S <= 1` trim
reduces to the Stage-A form. Part 2 drove the flow-coupled paths at `max|u| = 0.2009` (against
exactly `0` in part 1) and they stayed no-ops: the worst residual, `~9e-9` in ice volume at `t=10`,
is round-off amplified by convection, not a model difference — the three signatures below.

**Verdict: the Stage-A validation carries over to the Stage-B build.** B.1–B.5 are cleared to
proceed on `parties.stageB`. Artifacts: `ECCO_TESTS/StageB1_development/StageB_noninvasive/`.

**Part 2 — and why the residual is roundoff, not physics.** Compared only at the five times where
*both* runs have real data (the Stage-A reference writes every 2.0 time units; interpolating it
onto intermediate times manufactures an apparent 6.5e-3 "error" that is purely interpolation of a
curved `V(t)` — an early version of this analysis fell into exactly that trap):

| `t` | rel. diff ice volume | max pointwise `|ΔC_L|` | max pointwise `|Δθ|` |
|---|---:|---:|---:|
| 2 | 1.1e-12 | 3.5e-9 | 1.3e-8 |
| 6 | 7.0e-11 | 6.6e-9 | 1.9e-7 |
| 10 | 9.2e-9 | 6.4e-7 | 3.7e-5 |

Three independent signatures identify this as roundoff amplified by the convective flow's own
chaos rather than a model difference:

1. **Sign-random, not biased.** `mean(ΔC_L)/mean|ΔC_L|` is 0.03, 0.16, 0.37, 0.02, 0.23 across the
   five times — small and non-monotonic. A systematic effect (the IBM trim shaving ice, say) would
   drive that ratio toward 1 and hold it there.
2. **Exponential growth from a roundoff seed.** The difference grows at **0.652 per time unit**
   (doubling time 1.06) starting from `mean|ΔC_L| = 3.6e-12` at `t = 2` — a Lyapunov rate, not an
   offset. This is expected: the freshwater 288² case is the *vigorously convective* one, which is
   the entire basis of the Stage-A resolution argument (§2 of the validation summary). Differing
   struct layout and instruction scheduling between the two binaries perturb the last bit, and
   convection amplifies it.
3. **Against the project's own established bar it is negligible.** The accepted Stage-A binary
   lineage test (`aac878bf ≡ 6a8fa6eb`) recorded 1.7e-12 at `t=2` and **5.8e-5 at `t=6`** against a
   **7.2e-3** discriminating size. The Stage-B value at `t=6` is **7.0e-11** — roughly **six orders
   of magnitude tighter than the already-accepted lineage test** at the same time, and eight below
   the discriminating size.

**Stated precisely:** the Stage-B build is not bit-identical to Stage A, and cannot be — the build
is not byte-reproducible. What is established is that it is **indistinguishable from a rebuild of
the same physics**, by the same measure and against the same threshold the Stage-A campaign used to
certify its own binary lineage. The enthalpy drift is identical to four digits (1.201e-6 in both).

**Operational finding — every Stage-B deck must ship `p_fixed.inp` and `p_mobile.inp`.** With
`LAG_PARTICLE_RESOLVED` compiled in, PARTIES aborts at startup with
`Could not open p_fixed.inp` (`ParticleInput.c:87`, via `Particle_initialize`, `Cart3d.c:532`) even
when there are zero particles; both files need only contain `0`. This killed the first attempt at
part 2. It joins `stop.inp` (absent ⇒ abort at iteration 2, `Temporal_int.c:495`) and the module
set (`gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8` + `srun`) on the list of things that silently kill a job
at startup.

### B.0.1 Where runs must write — **`/anvil/scratch/x-mjalabert`, never `$HOME`**

`$HOME` is quota-limited at **25 GB** and was at **12.2 GB (48.9 %)** on 2026-08-09,
with **3.8 GB of that ECCO run output** because the early Stage-B jobscripts wrote
their run directories under `SLURM_SUBMIT_DIR`. Scratch is 511 GB of a 100 TB
allocation (**0.5 %**) — effectively unlimited.

The Stage-A campaign already did this correctly (`/anvil/scratch/x-mjalabert/YangFresh_07272026`
etc.); the Stage-B decks regressed. **Every Stage-B jobscript must put its run
directory on scratch** and keep only decks, analysers and summaries in `$HOME`:

```bash
RUN=/anvil/scratch/x-mjalabert/ECCO_StageB/<case>
mkdir -p "$RUN"
cp "$SLURM_SUBMIT_DIR"/{parties.inp,stop.inp,p_fixed.inp,p_mobile.inp} "$RUN/"
cp /home/x-mjalabert/PARTIES/PARTIES/parties.<binary> "$RUN/parties"
cd "$RUN" && srun ./parties > run.log 2>&1
```

Note `$HOME` is also carrying **4.2 GB of IDE server caches**
(`.vscode-server` 3.0 GB, `.antigravity-server` 1.2 GB) — the largest single item
after PARTIES itself, and safe to clear if space gets tight.

### B.1.0 NON-INVASIVENESS CONTRACT — every ECCO change is flag-gated

**Requirement (user, 2026-08-09): the ECCO work must not affect the code as it
stood before this project.** Audited 2026-08-09; this table is the contract and
must be re-checked whenever a new ECCO change lands.

| change | compile guard | runtime gate | reachable in a pre-ECCO run? |
|---|---|---|---|
| `Particle_release_by_interface` | `LAG_PARTICLE_RESOLVED && VOF_IBM` | `F_release <= 0` returns immediately | **no** — compiled out |
| `Particle_update_shell_liquid_fraction` | `LAG_PARTICLE_RESOLVED && VOF_IBM` | `F_release <= 0` returns | **no** |
| `Interpolate_integrate_shell_liquid_fraction` | `LAG_PARTICLE_RESOLVED && VOF_IBM` | — | **no** |
| motion lock + ramp (`Lagrangian.c:1133`, `:1467`) | `LAG_PARTICLE_RESOLVED && VOF_IBM` | `if (F_release > 0.0)` | **no** |
| `release.dat` diagnostic | inside `Particle_release_by_interface` | after the `F_release <= 0` return | **no** |
| `F_release` range check (`ParticleInput.c`) | `LAG_PARTICLE_RESOLVED && VOF_IBM` | only trips at `F_release >= 0.98` | **no** |
| `Conc_add_meltwater_RHS` | `PHASE_CHANGE && VOF_DIFFUSE` (**active in Stage A**) | call site `if (iconc == 2 && params->meltwater_tracer)`; default **0**; also needs `NConc >= 3` | **no** — runtime-gated off |
| scalar inits 34 / 35 | `CONC` (active) | new `switch` cases — an existing deck never selects them | **no** |
| VOF init 28 | none | pre-existing (Favier melting-RB), unchanged | n/a |
| sediment-clip restore (B.1.4-fix) | `VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE && VOF_IBM` | — | **no** — compiled out; and inert even if forced on, since `C_S ≡ 0` ⇒ nothing is ever trimmed in a sediment cell |
| `SED_CLIP_AUDIT` instrumentation | `DIFFUSE_TRACK_SEDIMENT_CLIP` (= `VOF_IBM && (RESTORE‖AUDIT)`) | — | **no** — compiled out, including its `MPI_Allreduce` |
| `diffuse_bound_liquid_fraction` split | none (**active in Stage A**) | second channel only computed when the caller passes non-`NULL` | **no** — every pre-existing call site goes through the `NULL` wrapper, and the fluid-side `MPI_Allreduce` is unchanged in count, type and op |

Stage-A flags are `#undef VOF_IBM` and `#undef LAG_PARTICLE_RESOLVED`, so the
entire release machinery is **compiled out**, not merely switched off. The only
additions on a Stage-A-reachable path are the meltwater tracer and the two new
init cases, and both are inert unless a deck explicitly asks for them.

**Defaults that must never be raised:**
`F_release = -1.0` (off) and `meltwater_tracer = 0` (off). Raising `F_release` to
the calibrated 0.709 as a *default* would silently start locking particles in
every existing deck — so the calibrated value lives in the ECCO decks and in the
`default.inp` comment block, **not** in the default itself.

Residual, stated honestly: new fields in `Parameters` and `Particle`
(`F_release`, `release_shell_cells`, `release_ramp_steps`, `meltwater_tracer`,
`phi_liq`, `t_released`, `release_ramp`) change those struct sizes even in a
Stage-A build. That cannot change behaviour — the fields are unread — but it does
mean the Stage-A binary is not byte-identical to the pre-ECCO one. The build was
never byte-reproducible anyway (§5 of the validation summary), so the contract is
enforced **functionally**, by the gate below.

**Verification: Stage-A no-impact gate ✅ PASSED 2026-08-09** (job 19761563,
binary `parties.stageA_recheck` `194a53b7`, built from the current source with the
ECCO flags OFF). Re-ran the recorded Le=100 salt gate:

| | `lambda` | rel. error | profile `L_inf` | salt drift | max speed |
|---|---:|---:|---:|---:|---:|
| recorded reference | 0.21844955 | 6.039e-03 | 1.441e-02 | 2.764e-03 | 0 |
| **current source, flags OFF** | **0.21844955** | **6.039e-03** | **1.441e-02** | **2.764e-03** | **0** |

Identical to every printed digit. **The ECCO work does not touch the pre-ECCO code
path.** Re-run this gate whenever a new ECCO change lands on a Stage-A-reachable
file (`Conc.c`, `VOF_DIFFUSE.c`, `Initial_Conditions.c`, `Velocity.c`).

### B.1 Additional code updates

#### B.1.1 Re-fit the EOS to polar conditions ✅ **DONE 2026-07-30**

Script: `ECCO_TESTS/StageB1_development/eos_polar_fit.py` (needs `gsw`; installed into the Anvil
anaconda python — note `pip` on the login node belongs to a *different* python 3.9,
so it must be invoked as `/apps/anvil/external/apps/anaconda/2025.06/bin/python3 -m pip`).

Least-squares fit of the **same** functional form the code implements
(`Velocity.c:2252`) to TEOS-10 over the polar box `SA ∈ [28,35] g/kg`,
`CT ∈ [−2.2, 2.0] °C`, `p = 0 dbar`, masked to points above the freezing line
(3152 of 3600). No code change — input values only.

| coefficient | Yang (lab) | **polar re-fit** |
|---|---:|---:|
| `Cb` [kg/m³/K²] | 0.011000 | **0.006935** |
| `T0` [°C] | 4.0000 | **3.6966** |
| `cS` [K/(g/kg)] | −0.2500 | **−0.2151** |
| `b0` [kg/m³/(g/kg)] | 0.7700 | **0.8116** |
| **RMS density error** | 1.10053 | **0.00143** |
| max abs. error | 1.39148 | **0.00452** |

A **770× reduction in RMS error**. At the box centre (`SA=31.5`, `CT=−0.5`) the
re-fit reproduces TEOS-10's expansion coefficients to `α: −0.66 %`,
`β: −0.01 %`.

**The fitted `Cb, T0, cS` are NOT physically meaningful individually.** The
density-maximum temperature `T_md(S) = T0 + cS·S` runs −2.33 °C (at SA=28) to
−3.83 °C (at SA=35), i.e. **below the freezing line everywhere in the box**
(−1.91 to −1.51 °C). The fit therefore samples only *one flank* of the parabola
and the three parameters are degenerate — Jacobian condition number **3.5e4**.
Only the local slope they imply is constrained. Quote `α, β`, not `Cb, T0, cS`.

**⚠ The scaling finding — this changes the ECCO deck, not just the EOS block.**
Carrying Yang's convention (`Δρ_ref = Cb·ΔT^q`, `betaT = 1`) into polar
conditions gives **`betaS ≈ 691`**. That is not a fitting artifact, it is the
physics: in polar water haline buoyancy dominates thermal buoyancy, the reverse
of Yang's warm freshwater tank.

| case | `α·ΔT` | `β·ΔS` (full meltwater) | ratio | `β·ΔS` (1 g/kg) | ratio |
|---|---:|---:|---:|---:|---:|
| ice shelf (ΔT=2.4 K, S=34) | 8.43e-5 | 2.67e-2 | **316×** | 7.84e-4 | **9.3×** |
| sea ice (ΔT=2.9 K, S=32) | 1.02e-4 | 2.51e-2 | **246×** | 7.84e-4 | **7.7×** |

With the thermal density scale the buoyancy field is `O(700)`, so the "free-fall"
velocity built from it is wrong by `√700 ≈ 26×` and `Ra_T` is **not** the
governing number. **ECCO must nondimensionalize on the haline scale**
`Δρ_ref = b0·S_m`, giving `betaS = 1` and `betaT = Cb·ΔT^q/(b0·S_m) = O(1e-3)`.
Same code, same equation, different input values — `q` and the `|·|` form are
untouched. `Re = √(Ra_S/Sc)` then, with `Ra_S = g(b0 S_m/ρ0)H³/(ν κ_S)`.

```
[eos]                        # ice shelf: T_i=-1.9 C, T_ocean=0.5 C, S_m=34
eos_q = 2.0
eos_betaT = 1.447514e-03
eos_betaS = 1.0
eos_Tmd0 = 2.331903
eos_Tmd_slope = -3.047358
```

**Validity range to state in the paper:** `SA ∈ [28,35] g/kg`,
`CT ∈ [−2.2,2.0] °C`, `p ≈ 0 dbar`; RMS 0.0014 kg/m³, max 0.0045 kg/m³; no
pressure/thermobaric terms. Thermobaricity is not negligible over a deep water
column (`α` rises **+4.3 %** at 50 dbar, **+17.2 %** at 200 dbar) but is
irrelevant across the cm-scale DNS box, which is the only vertical extent the
model resolves.

> **Open input:** `ΔT`, `S_m` and `H` above are placeholders for two plausible
> ECCO scenarios. The *dimensional* fit is scenario-independent and final; the
> nondimensional block must be regenerated once the ECCO case is fixed. Rerun
> `eos_polar_fit.py` with the case appended to `CASES`.

##### B.1.1a Independent cross-check, and one addition: prefer the LINEAR EOS (2026-07-30)

An independent re-derivation (`ECCO_TESTS/StageB1_development/eos_polar_refit.py`, written without reference to the
above) reproduces every qualitative conclusion of B.1.1 — the vertex degeneracy, the
below-freezing density maximum, and haline dominance (it gets 413× on its own box, against the
316×/246× tabulated above; the difference is only the choice of `ΔT` and `ΔS`). Two things it adds:

**1. The crossover is at `S ≈ 24 g/kg`, and Yang sits on the other side of it.**

| `S` (g/kg) | 0 | 10 | 20 | **24** | 28 | 34 |
|---|---:|---:|---:|---:|---:|---:|
| `T_freeze` (°C) | 0.00 | −0.52 | −1.06 | −1.28 | −1.50 | −1.85 |
| `T_maxdens` (°C) | 4.21 | 1.90 | −0.39 | **−1.31** | −2.21 | −3.55 |
| | inside | inside | inside | **crossover** | below | below |

Yang's `Sm = 5` lies well inside the range where the density maximum is *in* the liquid — which is
precisely the mechanism the benchmark's non-monotonic melt rate is about. ECCO lies beyond the
crossover, where meltwater is unambiguously buoyant.

> **Scope statement for the write-up:** the anomalous-convection **regime** validated in Stage A
> does *not* carry over to ECCO. What Stage A validated is the **machinery** — Stefan coupling,
> salt transport, Brinkman penalization, discrete enthalpy conservation. That is the correct and
> defensible claim, and it should be stated that way rather than implying regime transfer.

**2. A linear EOS is not merely adequate — it is the safer choice, because of extrapolation.**
With a free additive constant allowed (mandatory: buoyancy is a force, so a uniform offset is
absorbed into the pressure gradient — omitting it makes *every* model look ~99 % wrong as a pure
fitting artefact), over the ambient box:

| model | rms/std |
|---|---:|
| linear `c0 + c1·theta + c2·s` | 0.222 % |
| Roquet quadratic, vertex pinned to the TEOS-10 locus | 0.088 % |
| full quadratic (ceiling) | 0.001 % |

Both are far more accurate than needed *inside* the box. The decisive difference is **outside** it —
and the meltwater plume necessarily visits `S` far below 28 g/kg, since meltwater at the interface
is fresh. Error in the density anomaly against TEOS-10 (kg/m³), each model referenced to `S = 34`:

| `S` (g/kg) | 31.5 | 28 | 24 | 20 | 10 | 1 |
|---|---:|---:|---:|---:|---:|---:|
| quadratic | +0.003 | −0.000 | −0.015 | −0.040 | −0.148 | **−0.267** |
| linear | +0.006 | +0.009 | +0.008 | +0.002 | −0.033 | **−0.055** |

The quadratic wins inside the fit box and loses by **~5×** outside it, degrading monotonically as
the water freshens. Proportionately this is a "prefer linear" finding, not "the quadratic is
broken": 0.267 kg/m³ is ~1 % of the haline buoyancy scale `b0·ΔS ≈ 27 kg/m³`. But there is no
compensating benefit, since inside the box both are ~100× more accurate than required.

**Recommendation: for ECCO do not compile `EOS_NONLINEAR`; use the linear buoyancy path**
(`richardson`), with `d(b)/d(theta) = +6.45e-05` and `d(b)/d(s) = −2.66e-02` in the
`b = −(ρ−ρ_ref)/ρ_ref` convention — **verify the signs against the momentum assembly before
entering them**. If the quadratic is kept anyway, the B.1.1 block above is correct and safe
provided the fit box is extended to cover the dilution range actually visited.

#### B.1.2 Enable resolved IBM particle alongside VOF (`VOF_IBM`)
- **`Boundary.h`** — add `VOF_IBM`, `LAG_PARTICLE_RESOLVED`, `PARTICLE_RELEASE`.
- **`PARTIES/src/Eulerian/VOF_DIFFUSE.c`** (`VOF_DIFFUSE_compute_C_S`, ternary `C_L/C_S/C_G` Liu15
  mixing, ~line 1453) — reuse as-is; confirm no capillary/contact-angle path is active (`σ = 0`,
  `VOF_WETTING` off, MCL compiled out without `SURFACE_TENSION`).
- The melt-rate source (A.1.4) must skip the sediment region (`C_S`): no Stefan melting of the rock.
  **Already implemented** — `VOF_DIFFUSE.c:1344` skips any cell with
  `C_S >= DIFFUSE_SOLID_MASS_CUTOFF` (`= 0.05`, `VOF_DIFFUSE.c:21`). Without it the melt weight
  `w = C_L(1-C_L)` would be nonzero throughout the diffuse layer wrapping the *sediment*, turning
  the grain into a spurious melt/freeze source. Consequence to keep in mind for B.5: cells that are
  >=5 % sediment cannot melt, so a thin non-melting shell wraps the grain; at `d/Dx >= 24` this is a
  sub-percent shell, but it is the reason the release criterion is an average over the particle
  support rather than a single centroid sample.

#### B.1.3 Hold-in-ice → interface-triggered release

**Corrected design, 2026-07-30 — the particle starts in the FIXED list, not the release list.**

The original plan reused the existing `PARTICLE_RELEASE` machinery
(`Particle_release_to_mobile`, `Particle.c:864`, trigger `p->t_part_release < params->time`).
**That mechanism is wrong for ECCO.** `p_release_list` is *never painted*: every solid-indicator,
force-collection and MPI-update loop in `Lagrangian.c` iterates only `p_mobile_list` and
`p_fixed_list` (verified by inspection — `p_release_list` appears only in allocation
`Lagrangian.c:155`, initialization `Particle.c:202-229`, buffer sizing `Particle.c:264`, and the
transfer itself). A release-list particle therefore **does not exist in the flow until it
materializes**, which for a sediment grain locked in ice is physically wrong twice over: the grain
displaces no volume while frozen, and it appears from nothing at release.

The **fixed** list is the correct home. Fixed particles are painted into the solid indicator, and
`Lagrangian.c:762` computes their full hydrodynamic force set (`F`, `T`, `Int_U`, `Int_Omega`) every
step while their positions are *not* integrated. That is exactly "a rigid body present in the
domain, held stationary" — and it means that at the instant of release the particle already carries
a converged force history, so no cold-start shock arises (the problem the `STARTUP` block at
`Particle.c:401` exists to patch for time-based release).

**Second correction — no list transfer either. The grain is LOCKED IN PLACE, not moved.**

Transferring `fixed → mobile` mid-run is also wrong, for a mechanical reason: `p_list->Np` is the
**HDF5 dataset dimension** for particle output (`ParticleOutput.c:422`, `:704`). A transfer either
writes a short dataset (`Np` not updated) or changes the dataset dimension between snapshots (`Np`
updated) — both corrupt the output series.

The grain therefore **stays in the mobile list from `t = 0`** and its *motion* is locked. This is
physically identical — both a fixed particle and a mobile particle with `U = Omega = 0` paint the
same stationary rigid body into the flow — while leaving `Np`, the output layout and the MPI list
bookkeeping untouched. It also *improves* on the fixed-list route: because the grain is in the
mobile list throughout, it accumulates a converged hydrodynamic force history while still locked,
so at release there is no cold-start transient at all.

**✅ IMPLEMENTED 2026-07-30** — binary `parties.stageB_release` (`e0a8c8c1`, 1 086 464 B).

| piece | file | what it does |
|---|---|---|
| inputs | `IO/default.inp` | `F_release` (default **−1 = off**), `release_shell_cells` (2.0), `release_ramp_steps` (10) |
| params | `Include/DataTypes.h` | the three inputs; per-particle `phi_liq`, `t_released`, `release_ramp` |
| init | `IO/ParticleInput.c` | `t_released = −1` (locked) for **mobile** particles iff `F_release > 0`; fixed particles are never locked |
| shell integral | `Lagrangian/Interpolate.c` `Interpolate_integrate_shell_liquid_fraction` | local `∫(1−vf)·C_L dV` and `∫(1−vf) dV` over `r < R + release_shell_cells·h` |
| reduction | `Lagrangian/Particle.c` `Particle_update_shell_liquid_fraction` | one `MPI_Allreduce` of `2·Np` doubles → `phi_liq` on every copy |
| trigger | `Lagrangian/Particle.c` `Particle_release_by_interface` | frees a grain when `phi_liq >= F_release`; advances the ramp |
| lock + ramp | `Lagrangian/Lagrangian.c:1118` **and `:1455`** | zeroes `U, U_old, Omega, Omega_old` while locked; linear ramp after release — **both the predictor and the corrector halves**, see below |
| wiring | `Temporal_int.c` | trigger at the top of the step, shell integral beside `Interpolate_integrate_momentum`; **both gated on `which_stage == 0`** |

**Why the criterion is a shell average and not a centroid probe.** Inside the grain the ternary
mixture sets `C_S = 1` hence `C_L = 0`, so *any* sample of the particle interior reads "no liquid"
forever and the grain would never be released. The weight `(1 − vf_cell)` — the same level-set
volume fraction the IBM itself uses — is zero in cells the particle fills and one in cells it does
not, so the criterion measures exactly the melt state of the grain's *surroundings*. The finite
shell thickness matters too: the non-melting shell of §B.1.2 keeps ice slightly longer in the first
ring of cells, so a one-cell probe would read that artifact rather than the melt front.

**[MPI] Why this is safe by construction.** The shell straddles rank boundaries, so no single rank
can evaluate it. Each rank integrates only cells it owns; one `Allreduce` sums numerator and
denominator; every rank then divides *the same two numbers*. The decision is bit-identical on every
rank — there is no broadcast to get wrong and no way for one rank to free a grain another still
considers locked. It is also decomposition-independent (pointwise integrand, each cell visited
once), so it clears the rank-count test by design. Requires `LIST_STATE_BOTH` at the call site,
which is why the integral sits beside `Interpolate_integrate_momentum`; the degradation mode if a
rank held shell cells but no particle copy is graceful (numerator and denominator lose the same
cells, so the ratio stays a valid average over the cells visited).

**Semantics to keep in mind:** with `F_release > 0` *every* mobile particle starts locked and is
freed individually as its own surroundings melt. For ECCO (one grain) that is exactly right, and it
generalizes correctly to N grains. `F_release <= 0` disables the mechanism entirely and reproduces
previous behaviour exactly — which is what the default-off regression below tests.

**Interpretation of the threshold:** `F_release` is the *fraction of the grain's immediate
neighbourhood that is liquid*. A half-exposed grain is still mechanically keyed into the ice, so
values near 1 (0.9–0.95) correspond to "essentially free"; 0.5 releases a grain that is still
half-embedded. This wants calibrating in B.5 against when the grain actually loses contact.

**✅ Criterion verified independently of the solver (2026-07-30).** The shell integral was
replicated in numpy and evaluated on prescribed `C_L` fields, for the shakeout geometry
(`R = 0.04`, `d/Δx = 20.5`, shell = 2 cells, grain spanning `y = 0.62–0.70`):

| prescribed field | `phi_liq` |
|---|---:|
| grain fully in ice (`C_L = 0`) | **0.0000** |
| grain fully in water (`C_L = 1`) | **1.0000** |
| melt front at the grain's bottom, `y = 0.62` | 0.1378 |
| melt front at the grain's centre, `y = 0.66` | **0.4937** |
| melt front at the grain's top, `y = 0.70` | 0.8599 |
| melt front just above, `y = 0.71` | 1.0000 |

Correctly bounded in `[0,1]`, monotone in the front position, and `0.4937` at the half-way front —
i.e. it measures what it claims to. The discrete shell area matches the analytic annulus
`π(R_out² − R²)` to 1.8 %. With `F_release = 0.9` the grain is freed between the last two rows,
i.e. once the front has cleared its top by about one cell — the intended semantics.

**⚠ BUG FOUND BY THE SHAKEOUT — the lock must be applied TWICE (fixed 2026-07-30, job 19581346).**

The first 2-D shakeout showed the grain creeping downward at `U_y = −4.2e-4` while
`t_released = −1` (locked): `X_y` drifted `0.660000 → 0.659107` over `t = 0 → 3.5`. Small, steady,
and entirely wrong — a grain frozen in ice was settling.

Cause: **`Lagrangian_integrate_particle_motion` is a predictor–corrector.** The predictor loop
(`Lagrangian.c:1006`) computes `U` from the forces; then `Collision_evaluate` runs; then a **second
mobile loop** (`:1359`) *recomputes* `U` from `U_old` and the forces at `:1392` and integrates `X`
again. The lock was only in the predictor, so the corrector rebuilt a non-zero `U` from the
hydrodynamic force the grain still feels — and, because the predictor had zeroed `U_old`, it
rebuilt it from zero every step, giving exactly the small constant creep observed.

This was a **process failure, not just a coding slip, and it is worth recording as such.** When the
edit was first applied, the tool reported two identical candidate blocks and I assumed the second
belonged to the fixed-particle `STUP_INIT_PRIMPOSED` path — without checking. It was the mobile
corrector. The lesson: when a codebase presents two byte-identical blocks, identify both before
patching either.

A full audit of every `U`/`Omega` assignment in `Lagrangian.c` was then done to confirm there is no
third site: the remaining ones are the `OSCILLATING_PARTICLE` and `STARTUP` branches (both inside
the two loops, both upstream of the lock), the `LEFT_WALL_VELOCITY_FREESLIP` pin (not compiled), and
the `ubulk_target` blocks, which are in the **fixed**-particle loops and never touch a mobile grain.

**Restart persistence — a second bug, found by reading, before it cost a run.** The per-particle release
state is *not* in the particle HDF5 element list, and the MPI particle datatype is
`sizeof(Particle)` raw bytes (`Particle_initialize_MPI_datatype`, `Particle.c:290`) so rank
migration is safe — but **restart was not**. A resumed run would have found `t_released`
uninitialized and, on the locked branch, would have **re-locked a grain already in free fall and
re-ramped its velocity from zero**. That is a physics error, not a lost diagnostic, and it would
have been easy to mistake for a collision or a drag artefact.

Fixed on three sides:
- `ParticleOutput.c` now writes `t_released` (one more `OUTPUT_ELEMENT`, plus its case in the
  element switch).
- `ParticleInput.c` reads it back **guarded by `H5Lexists`**, so restart files written before
  B.1.3 still resume; when the field is absent it says so on rank 0.
- `Particle_initialize_nonessential_data` defaults every grain to **locked** (`t_released = -1`).
  Locked is the conservative default: a grain wrongly held is visible immediately, whereas a grain
  wrongly freed silently corrupts the trajectory.

**✅ Default-off regression PASSED (job 19581196, 2026-07-30).** The Le=100 gate re-run on
`parties.stageB_release` with `F_release` at its default:

| | λ | rel. error | profile L∞ | salt drift | max speed |
|---|---:|---:|---:|---:|---:|
| Stage-A reference | 0.21844955 | 6.039e-03 | 1.441e-02 | 2.764e-03 | 0 |
| Stage-B flags only (`85e0d33f`) | 0.21844955 | 6.039e-03 | 1.441e-02 | 2.764e-03 | 0 |
| **+ B.1.3 release code (`e0a8c8c1`)** | **0.21844955** | **6.039e-03** | **1.441e-02** | **2.764e-03** | **0** |

Identical to every printed digit across all three. The release machinery is a verified no-op when
disabled, so it cannot have perturbed anything validated in Stage A.

#### B.1.3-fn Does the release actually FIRE?  — **two MPI bugs found and fixed 2026-07-30**

The default-off regression above proves the release code is a **no-op when
disabled**. It says nothing about whether it *works when enabled*, and until that
is shown the mechanism is not ECCO-relevant. Deck: `ECCO_TESTS/StageB1_development/StageB_release_fn/`
(ice above, ocean below, warm bottom wall, one grain `d/Δx = 25.6` frozen half a
diameter above the melt front; analyser `analyze_release.py`, which reads its
constants from the deck so it runs against any release case).

**A `release.dat` diagnostic was needed first.** The only observable was a
one-shot stdout line, which cannot calibrate `F_release` — that needs `φ_liq(t)`
*through* the crossing, plus `U(t)` across it to show the handoff is impulse-free.
The new file is rank-0, Stage-B-only, appends on resume, and changes no existing
output format.

Building it immediately exposed two defects, **both invisible in the default-off
regression** because that test has no particle at all:

| # | defect | consequence |
|---|---|---|
| 1 | `release.dat` written by walking rank 0's **local** particle list | **Rank 0 usually holds no copy of the grain** (the domain is split; the grain sits wherever it sits), so its list is empty and the file got a header and *zero data rows*. Confirmed empirically: 54 bytes, header only. |
| 2 | the `[release]` stdout notice sits **inside** that same loop under `rank == 0` | **Pre-existing bug in B.1.3 itself, not in the diagnostic.** On any decomposition where rank 0 does not own the grain, the release fires and leaves **no trace in the log whatsoever**. |

Both are fixed by reducing over ranks before writing: each rank contributes its
copies to an `MPI_Reduce`, and rank 0 divides by the copy count. Copies carry
identical values (`φ_liq` comes from the existing `Allreduce`; `t_released` and the
ramp are decided identically on every rank from that same number), and a sum of
zeros over any count is still exactly zero — which is what the "locked" check reads
as exact. The announcement is emitted from the reduced values.

Binary `parties.stageB_rel3` (`9908667b`). Run 19581597.

> **Note on the parallel `StageB_release2d` run (job 19581346).** It uses
> `parties.stageB_full`, built before these fixes, so it will produce a header-only
> `release.dat` and print no release notice. Its **physics is unaffected** — both
> bugs are in reporting only — and the release instant is still recoverable
> post-hoc from `mobile.dat`, where `U` departs zero. It is not a wasted run.

#### B.1.3-fn RESULT — ✅ **PASS 2026-07-30** (job 19581767, `parties.stageB_rel4` `993a2972`)

After the three fixes, all four checks pass on the 2-D freshwater deck
(`ECCO_TESTS/StageB1_development/StageB_release_fn/run32b`, 25 018 steps to `t = 25`):

| check | result |
|---|---|
| release fires | ✅ `t = 19.026013`, `φ_liq = 0.9000` at threshold `0.9` |
| **locked before release** | ✅ `max|U| = 0.000e+00`, `|Δy| = 0.000e+00` — **exactly** zero for 19 026 steps |
| impulse-free handoff | ✅ peak `|dU_y/dt| = 0.2593` = **0.52×** free-settling `a_free = 0.5` (bar: <2×) |
| settles after release | ✅ falls, mean `U_y = −1.0e-3` |

**The handoff is a spin-up, not a kick.** `U_y` is exactly 0 through the lock and
the first two ramp steps, then grows monotonically over the 10 ramp steps
(`−2.0e-5 → −5.6e-4`); the peak acceleration lands **one step after the ramp
closes** and decays as drag balances. An unramped release would have put the peak
*at* the release instant.

**The third bug — the lock leaked, and it was in the physics.**
`Lagrangian_integrate_particle_motion` is predictor/corrector. The lock was applied
to the predictor only, so the corrector rebuilt a non-zero `U` from the force the
grain still felt and it **crept downward while nominally frozen in the ice**:
measured `max|U_y| = 1.44e-4` and `−3.7e-5` of drift (0.0095 cells) in `t = 0.44`.
Fixed by applying the same lock and ramp in both halves (`Lagrangian.c:1133` and
`:1467`). The other two position updates (`:1306`, `:1599`) are `STARTUP`-gated
fixed-list paths and are not on this path.

> **Why the earlier evidence was insufficient.** Three defects, none reachable by the
> default-off regression, which has no particle at all: two made the release
> **unobservable** (empty `release.dat`; the `[release]` notice never printing), and
> the third made it **wrong** (pre-release creep — a grain that drifts before it is
> freed has the wrong release position and a spurious pre-release wake). "No-op when
> disabled" was never evidence that the mechanism works when enabled.

##### `F_release` calibration — answers the question B.5 deferred

`φ_liq` measured against how far the melt front has passed the grain (0 = front at
the grain's underside, 1 = grain just fully uncovered):

| exposure | 0.5 (half uncovered) | **1.0 (just uncovered)** | 1.5 | 1.67 |
|---|---:|---:|---:|---:|
| `φ_liq` | 0.565 | **0.709** | 0.867 | 0.900 |

**`F_release = 0.9` fires 0.67 diameters LATE** — the front has already cleared the
whole grain. For "released the moment it is mechanically free", use
**`F_release ≈ 0.71`**.

> **⚠ Hard upper bound: `φ_liq` saturates at ≈ 0.977, never 1.** The shell always
> contains part of the diffuse CH band and the non-melting sediment shell of §B.1.2.
> **`F_release ≥ 0.98` never releases the grain at all.** This is a silent-hang trap
> and must be range-checked.

#### B.1.4 RESULT — tracer is CORRECT; the 9.3 % is an **uncompensated clip in the IBM coupling**

**Discriminator run (job 19761565, `StageB_tracer_noparticle`): the parallel
`StageB_release2d` deck with the particle removed and nothing else changed.**

| `t` | no particle | with particle |
|---|---:|---:|
| 2 | −1.08e-7 | −7.20e-2 |
| 10 | −4.77e-7 | −6.75e-3 |
| 15 | **−5.35e-7** | +2.29e-2 |
| 30 | — | +9.30e-2 |

**Without the particle the identity `∫C_mw = ∫dF` holds to 5.4e-7** — five orders of
magnitude inside the 1e-3 bar. So the tracer's transport, its source term and the
RK time integration are all correct.

> **My earlier hypothesis was wrong.** I proposed an RK quadrature mismatch between
> the F equation's `ch_rhs_n/nm1` history and the tracer's `melt_src/melt_src_old`.
> The no-particle run refutes it: with the identical integration and no particle,
> the drift vanishes. The defect is entirely IBM-coupled.

**Mechanism, located in the code.** `diffuse_bound_liquid_fraction`
(`VOF_DIFFUSE.c:1110`) trims `C_L → 1 − C_S` in **every** cell, but accumulates the
removed mass for redistribution **only** where `cs < DIFFUSE_SOLID_MASS_CUTOFF`
(`:1127`):

```c
vof->C_L[k][j][i] = new_cl;                     /* trim applied everywhere   */
if (cs < DIFFUSE_SOLID_MASS_CUTOFF) {           /* ...but only counted here  */
    local_delta += (new_cl - old_cl) * measure;
}
```

The caller then returns `local_delta` to the interfacial band, which is what makes
the clip conservative — and the surrounding comment says exactly why that matters
("a plain clip is a hidden enthalpy leak that biases the Stefan front speed,
measured ~5 % in the A.3 gate"). **Liquid trimmed inside sediment cells is exempt
from that accounting and is silently destroyed.**

Quantitatively consistent with a continuous per-step leak, not a startup artifact:
at `t = 30` the deficit is **462 cell-units against a particle volume of 171**
(2.71×), through **332 cells** with `C_S ≥ 0.05`, while `C_S` itself is *exactly*
unchanged (`max|ΔC_S| = 0`) because the grain is locked and sediment cannot melt.
The early **negative** drift is the `t = 0` transient: init 28 sets `C_L = f`
everywhere including inside the grain, and the first trim clips that away.

**Consequences.**
1. The meltwater tracer is **validated** for entrainment work in the absence of a
   resolved grain, to 5e-7.
2. `∫dF` is **not** a valid proxy for `∫m` once a resolved particle is present, so
   the "9.3 % tracer failure" was the *check* being invalid, not the tracer.
3. **But the leak is real and is not a diagnostic artifact** — it is genuine liquid
   mass destroyed in the ice/sediment overlap, and it biases the ice-volume budget
   that every ECCO melt number rests on. It grows with the sediment surface area,
   so 3-D and multi-grain cases will be worse.

#### B.1.4-fix DECISION AND IMPLEMENTATION — option (a), rim-local — 2026-08-28

The three options were (a) redistribute the trimmed liquid into the band like all
other clipped mass, (b) prevent it entering the solid by tightening the IBM forcing,
(c) leave it and report a known budget error.

**Chosen: (a), refined to deposit rim-locally rather than globally.** Reasons, in
the order they decide it:

- **(c) is not available.** The leak biases the ice-volume budget every ECCO melt
  number rests on, and it scales with sediment surface area — so the configuration
  where it is worst is exactly the production configuration (3-D, and later
  multi-grain). A known-error caveat on the central quantity is not a deliverable.
- **(b) is misdiagnosed.** The IBM constrains *velocity*; the transport that carries
  `C_L` into the grain is the Cahn–Hilliard mobility flux in the implicit biharmonic
  solve, which the IBM never touches. Tightening the forcing cannot close it. The
  operator-level version of (b) — masking the CH mobility by `C_S`, the way the salt
  diffusivity is masked by `F` in G.5 — *would* be correct, but it converts the
  constant-coefficient biharmonic into a variable-coefficient one, changing a
  validated Stage-A implicit operator and its PCG conditioning, and would need its
  own gate. Wrong risk to take at this point; recorded as the fallback if (a) proves
  insufficient.
- **(a) is what the code's own mass definition already implies.**
  `diffuse_liquid_mass_integral` (`VOF_DIFFUSE.c:1092`) *already* excludes sediment
  cells from `∫C_L`. Under that definition liquid diffusing into a sediment cell has
  left the accounting domain, so `∫C_L` over the fluid region drops with no physical
  cause unless the mass is returned. Returning it is bookkeeping consistency, not a
  new modelling assumption.

**The refinement, and why it matters.** Naive (a) hands the mass to
`diffuse_redistribute_liquid_mass_delta`, which is *globally* collective — it
redeposits anywhere on the interface with `F(1−F)` headroom. The `VOF_DIFFUSE_CLIP_AUDIT`
comment block already warns what that costs in 2-D (row-melt CV 0.35–0.64: teleporting
melt between scallop rows flattens the interfacial dynamics that set convective heat
delivery). So the sediment channel deposits into the **diffuse rim of `C_S` outside
the grain surface** — cells with `0 < C_S < DIFFUSE_SOLID_MASS_CUTOFF` — which is
precisely where the spurious flux took the mass from. The global routine is kept only
as a fallback when that rim has no headroom left, and the residual it has to absorb is
reported, so "how much had to be teleported" is a measured number rather than an
assumption.

The rim is well-defined and several cells thick by construction:
`VOF_DIFFUSE_compute_C_S` builds `C_S = ½ − ½tanh((r − (R − shift))/denom)` with
`shift = √2·ln(19)·Cn` and `denom = 2√2·Cn`, which places `C_S = 0.05` *exactly* on
`r = R`. So `0 < C_S < 0.05` is the shell `r > R`, out to the solid support.

**Implementation (`VOF_DIFFUSE.c`).**

| piece | what it does |
|---|---|
| `diffuse_bound_liquid_fraction_split(db, &sed)` | the clip, with the trim accounted in two channels: fluid cells (return value, the Stage-A quantity) and sediment cells (`*sed`, previously discarded). `diffuse_bound_liquid_fraction(db)` is now a wrapper passing `NULL`, so all five existing call sites are untouched. |
| `diffuse_sediment_rim_eligible(cl, cs, pass)` | host-cell predicate: `0 < C_S < cutoff`, two passes (strict then relaxed band), mirroring `diffuse_mass_redist_eligible`. |
| `diffuse_redistribute_sediment_clip_delta(db, delta)` | rim-local deposit, `F(1−F)`-weighted, headroom-capped; returns what it could not place. |
| wiring | both in-loop clip sites: the main post-biharmonic clip, and the `Pi_adm` trim after the contact-angle projection (a **second, previously unnoticed leak channel** — that call site also discarded its return value). |
| `VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT N` | the instrumentation the roadmap asked for: `SED_CLIP_AUDIT` lines every `N` steps with `sed`, `sed_resid`, their running sums and `cum|sed|/cum_melt`. Log-only. |

**Non-invasiveness — by construction, not by argument.**
`DIFFUSE_TRACK_SEDIMENT_CLIP` is `1` only when `VOF_IBM && (RESTORE || AUDIT)`, so a
Stage-A build compiles out the entire sediment channel including its `MPI_Allreduce`.
The fluid-side reduction is deliberately left as a separate 1-element `MPI_Allreduce`
rather than folded into a 2-element one: MPI reductions are deterministic per
`(count, datatype, op, communicator)` but **not across counts**, and the Stage-A gate
has to reproduce `λ` to every printed digit — not worth spending that guarantee to
save one collective on the ECCO path. Both new flags default `#undef` in `Boundary.h`.

**The `t = 0` transient is correctly excluded.** Init 28 paints `C_L = f` inside the
grain and the first trim clips it away — mass that is *not* physical and must not be
redistributed. It is excluded automatically: the two initialization-time
`diffuse_bound_liquid_fraction` calls (`VOF_DIFFUSE.c:576`, `:837`) run after `C_S` is
built and discard their return, so only clips inside the time loop are accounted.

**[MPI]** pointwise integrand, one `Allreduce` for the weight and one for what was
applied; `delta` is derived from reduced values only, so every rank runs identical
loop trip counts and the collectives cannot desynchronize. Decomposition-independent
by the same argument as B.1.3's shell integral.

**Verification submitted 2026-08-28** (see B.1.4-fix RESULT below when it lands):

| job | arm | purpose |
|---|---|---|
| 20198502 | `parties.stageA_b14`, ECCO flags OFF | B.1.0-mandated Stage-A no-impact gate — B.1.4 touches `VOF_DIFFUSE.c`, so this is required. Must reproduce `λ = 0.21844955`. |
| 20198503 | `parties.stageB_sedctl` — audit ON, restore **OFF** | control: measures the leak directly instead of inferring it from the tracer |
| 20198504 | `parties.stageB_sedfix` — audit ON, restore **ON** | identical deck, one compile flag apart |
| 20198505 | `parties.stageB_sed3d` | B.2 3-D shakeout, fix ON (see §B.2) |

Build script: `ECCO_TESTS/StageB1_development/build_b14.sh` (patches `Boundary.h`, builds all four,
restores `Boundary.h` under a `trap`). Decks: `ECCO_TESTS/StageB1_development/B14_sedclip/{ctl,fix}`.

**Pass criteria for the A/B.** (1) `cum|sed_resid| ≈ 0` — the rim absorbed the mass
locally and nothing had to be teleported; (2) the tracer identity `∫C_mw = ∫dF` closes
to ≤ 1e-3, against the control reproducing the recorded +9.3e-2 at `t = 30`;
(3) release time and settling trajectory not materially disturbed.

#### B.1.4-fix ROUND 1 RESULT — ❌ **THE FIRST IMPLEMENTATION WAS WRONG. The A/B caught it.** (2026-08-29)

The **decision** (option (a) — restore rather than discard) survives unchanged. The
**host set** chosen for the deposit did not, and the failure was severe.

| job | result |
|---|---|
| 20198502 Stage-A gate | ✅ **PASS** — `λ = 0.2184495545645199`, bit-identical to the recorded value in all 16 digits, as are `fit_rms`, `profile_L∞`, salt drift `2.7638540325809822e-03`, enthalpy `4.975384548799866e-10`, `max_speed = 0`. The clip-accounting split and the non-invasiveness machinery are sound. |
| 20198503 control | ✅ **faithful** — `max φ_liq = 0.66812110287866755`, identical to **17 digits** to the recorded Jul-30 baseline. The instrumentation is genuinely log-only. |
| 20198504 fix | ❌ **destroyed the ice slab.** |
| 20198505 B.2 3-D | ⚠️ **contaminated** (ran the broken binary); see §B.2. |

**What the fix arm did.** Far-field ice `C_L` at `y ∈ [0.80, 0.85]`, nowhere near the
grain:

| | `t = 2` | `t = 10` | `t = 30` |
|---|---:|---:|---:|
| control | 0.0000 | 0.0000 | 0.0008 |
| **fix** | 0.0057 | **0.9539** | **1.0000** |

The entire ice slab melted. Consequently the grain's shell liquid fraction ran to
0.976 and it released at `t = 1.98` — the control (and the recorded baseline) never
reach the 0.9 threshold at all within `t = 30`.

**Root cause — a locality argument keyed on a tanh tail.** §B.1.4-fix argued the host
set `0 < C_S < DIFFUSE_SOLID_MASS_CUTOFF` "is the shell `r > R`, a few cells thick".
The first half is right: `VOF_DIFFUSE_compute_C_S`'s `shift` does place `C_S = 0.05`
exactly on `r = R`. **The second half was never checked and is false.** `C_S` is a
tanh whose support runs to `diffuse_solid_support_extra(Cn) = −shift + denom·atanh(1 − 2·10⁻¹²)`,
which for this deck is **0.102 against `R = 0.04`**. So the "rim" was a disc of radius
`3.5 R` — **3822 cells, 5.8 % of the domain**, against the ~180 a true 2-cell shell
would hold. `F(1−F)` weighting over that region, renormalised by a tiny total weight,
then concentrated the delta into whatever few cells passed the band test and saturated
them; repeated every RK stage for 12 000 steps it converted the ice to liquid.

**A second defect, in the instrumentation.** The residual was reported *after* the
global-fallback call, so it could not distinguish "placed locally" from "teleported
across the domain". This is exactly why the 3-D run looked healthy (`sed_resid ~ 10⁻²⁰`,
`ratio = 1.4 × 10⁻²`) while doing the wrong thing: with the grain buried in ice the
adjacent cells held no liquid, the local pass placed nothing, and the global routine
quietly smeared the mass over the whole melt interface.

#### B.1.4-fix ROUND 2 — corrected design (submitted 2026-08-29)

| change | why |
|---|---|
| host set = **one-cell dilation of the grain** (`diffuse_sediment_adjacent`): non-sediment cells with a face neighbour at `C_S ≥ cutoff` | locality now keyed on the **mesh**, not on a tanh level set. Needs no halo exchange: `compute_C_S` fills the ghosted range `L_*` from particle geometry on every rank intersecting the support, so `C_S` ghosts are already valid; only owned cells are written |
| weight by **`C_L`**, and skip cells with `C_L ≤ 0.005` | the flux can only have come from neighbours that *had* liquid; a bulk-ice neighbour must not receive any |
| **refuse** the delta when no adjacent cell holds liquid | if nothing next to the grain is liquid there was no leak to undo — a non-zero trim then is `t=0` transient or round-off, and amplifying it is the failure mode above |
| per-cell per-stage ceiling `DIFFUSE_SED_MAX_DCL = 0.05` | a one-face, one-stage diffusive leak can never legitimately owe a cell an O(1) change in `C_L`. Kills the concentration mechanism outright |
| **global fallback removed** at both call sites; unplaceable mass is **dropped and reported** | teleporting is the hazard `CLIP_AUDIT` warns about, and the silent fallback is what made a broken restore look healthy in 3-D |
| audit reports **signed net** and **dropped** separately | `Σ|trim|` overstates the destruction badly (the per-stage trims largely cancel: control `cum\|sed\|/cum_melt = 0.44` against a *net* tracer drift of 9.3 %). Net is the physical number; a large `dropped` is the signal to escalate |

**Escalation path if `dropped` is not small:** mask the CH mobility by `C_S`, the way
G.5 masks salt diffusivity by `F`. That is the operator-level correct fix; it was
deferred because it converts the validated constant-coefficient biharmonic into a
variable-coefficient one and needs its own gate.

**Round-2 jobs: 20213564 (Stage-A gate), 20213565 (control), 20213566 (fix).**
The 3-D shakeout is **deliberately not resubmitted concurrently this time** — running
it against an unvalidated fix already cost 12 h × 64 ranks. It is gated on this A/B.

**Round-2 Stage-A gate ✅ PASSED** (20213564): `λ = 0.2184495545645199`, again identical
to the recorded value in all 16 digits. The round-2 edits are confined to the ECCO path.

> **Operational lesson — every jobscript now gets a job-id-unique scratch directory.**
> The B.0.1 jobscript template opens with `rm -rf "$d"; mkdir -p "$d"` on a fixed path,
> so resubmitting an arm **destroys the run it is meant to be compared against**. That
> is what happened here: the round-2 submissions wiped the round-1 output while it was
> still being analysed, and the round-1 *raw fields* are unrecoverable. The round-1
> findings survive only because the decisive numbers had already been extracted into
> the table above — luck, not process. All four ECCO jobscripts now use
> `.../<case>_${SLURM_JOB_ID}` and `mkdir -p` without `rm -rf`. Anything comparative
> should follow that pattern.

#### B.1.4-fix ROUND 2 RESULT — ⛔ **option (a) is empirically dead in BOTH its forms** (2026-08-29)

Jobs 20213565 (control) / 20213566 (fix), both to `t = 30`.

| | control (restore OFF) | fix (rim-local restore ON) |
|---|---:|---:|
| far-field ice `C_L`, `t=30` | 0.0008 | 0.0007 |
| ice volume, `t=30` | 0.32680 | 0.32138 |
| **ring `C_L`** (`R < r < R+2Δx`), `t = 2 / 10 / 30` | 0.258 / 0.416 / **0.675** | 0.507 / 0.967 / **0.993** |
| `max φ_liq` | **0.6681** | **0.9769** |
| release | never fires | fires at **`t = 9.07`** |
| audit `net` / `cum_melt` | −2.138e-4 / 4.854e-4 = **44 %** | −2.018e-4 / 4.697e-4, `dropped` 1.041e-4 (**52 %**) |

**Round 2 fixed the catastrophe and exposed the real problem.** The ice slab now
survives (round 1 melted it entirely), the mesh-based host set and the per-cell ceiling
work as intended, and the Stage-A gate is still bit-identical. But the restore pumps
liquid into the ring around the grain until it **saturates at `C_L ≈ 0.99`** — a water
bubble around a rock buried in ice, with no heat source to have melted it — and that
false wetting trips the release criterion at `t = 9.07`.

**Neither arm is right, and the audit's `net` column says why.**
`net ≈ −cum|sed|` to four digits in both arms: **every trim is a removal, none cancel.**
(The comment written in round 2 claiming the per-stage trims "largely cancel" is wrong
and is corrected here.) So `C_L` is being regenerated inside the grain *every stage* at
~2.3e-4 per grain cell, drawn by the `C_L ≈ 0.68` ring across one cell of the
Cahn–Hilliard operator. That is a genuine diffusive flux, and it leaves only two
outcomes:

- **destroy it (control)** ⇒ the grain is a permanent *sink*: it eats liquid from its
  surroundings for the whole run, 44 % of the melt budget. The ring is left artificially
  **dry**, which is why `φ_liq` stalls at 0.668 and the grain **never releases** — in
  2-D here, and in the 3-D shakeout out to `t = 125`.
- **give it back (fix)** ⇒ conservation improves, but the returned mass has nowhere
  physical to go, so the ring saturates and the grain releases **spuriously early**.

**Conclusion: the clip is a symptom; the disease is that the CH operator does not know
the sediment is there.** Any post-hoc bookkeeping on the trimmed mass produces an
artifact, because the flux should never have existed. This is now an *empirical* result,
not the argument-from-first-principles that B.1.4-fix used to defer the operator-level
option — and it settles that question against option (a).

**This also blocks B.2 either way**, which is the operationally important part: the
control suppresses the release (nothing to shake out), and the fix triggers it
spuriously (meaningless to shake out). §B.2 cannot produce a valid result on either
build.

#### B.1.4-fix ROUND 3 — CH mask, Dirichlet form — ✅ **killed the leak, ❌ wrong boundary condition** (2026-08-29)

Decision taken with the project lead: do **both** the operator-level fix and a
re-basing of the release criterion, in that order.

`VOF_DIFFUSE_SEDIMENT_CH_MASK` gives sediment cells an **identity row** in the
Cahn–Hilliard biharmonic solve, with `C_L` and the RHS zeroed there once before the CG
(every iterate then stays zero, so no per-iteration masking is needed and the stencil
reads the correct value). Symmetric, so PCG is unchanged; inert without `VOF_IBM`.
`VOF_DIFFUSE_SHELL_FRACTION_NONSOLID` renormalises the B.1.3 criterion by the non-solid
shell volume, `φ = ∫(1−vf)C_L / ∫(1−vf)(1−C_S)`.

| job | | result |
|---|---|---|
| 20215987 | Stage-A gate | ✅ `λ = 0.2184495545645199`, enthalpy `4.975384548799866e-10` — bit-identical again |
| 20215988 | mask, old shell norm | leak **eliminated**: `net = 0.000000e+00`, `cum\|sed\| = 0`. `cum_melt = 4.859e-4` vs control `4.854e-4` (0.1 %). But `max φ_liq = 0.082` |
| 20215989 | mask + non-solid norm | same, `max φ_liq = 0.084` |

**The mask does exactly what it was built to do** — no liquid is ever transported into
the grain, so there is nothing to trim, nothing to restore, and neither round-1 nor
round-2 artifact can exist. The melt budget is untouched.

**But the Dirichlet form is the wrong boundary condition, and it shows.** `C_L` against
distance from the grain surface at `t = 30`:

| `r − R` (cells) | −1→0 | 0→1 | 1→2 | 2→4 | 4→8 | 8→16 |
|---|---:|---:|---:|---:|---:|---:|
| control | 0.640 | 0.680 | 0.671 | 0.692 | 0.672 | 0.644 |
| **mask (Dirichlet)** | 0.040 | 0.033 | 0.100 | 0.255 | 0.440 | **0.536** |

A **depletion layer ten cells deep**. `C_L = 0` in the rock asserts that the rock forces
*ice* against itself, which is not a modelling choice anyone made — it is an artifact of
using a zero-*value* condition where a zero-*flux* condition belongs. The shell
renormalisation could not rescue it (0.082 → 0.084) because by then `C_L` genuinely is
~0 near the grain; the normalisation was never the binding constraint here.

#### B.1.4-fix ROUND 4 — zero-FLUX at the grain surface (submitted 2026-08-29)

The rock is **indifferent** between water and ice: neither phase crosses into it, and
neither is preferred. That is homogeneous **Neumann** on `C_L` at the sediment boundary
— equivalently the 90° contact angle the deck already specifies (`contact_angle_deg = 90`)
and consistent with B.1.2's "no capillary/contact-angle path active".

`TwodOps_scalar_laplacian` is shared with the pressure solver and must not be touched, so
rather than rewriting the stencil, `diffuse_sediment_noflux_correction` **cancels the
face terms that reach into sediment**. Each face enters the 5/7-point stencil as
`(φ_nb − φ_c)·idx_c·idx_u`, so subtracting exactly that term makes the face contribute
zero — an exact no-flux face, not an approximation. The identity rows are kept, so the
grain is a **hole in the CH domain with no-flux walls**. Symmetry survives (the fluid row
drops its coupling to a sediment cell that never coupled back), and `L∘L` composes two
such symmetric operators, so PCG remains valid. Guarded `#error` against `AXISYM_RZ`,
which has extra `r_u` metric factors this correction does not yet carry.

> **Deck correction — "does the release fire" was never testable on this deck.** The
> control's `φ_liq` at `t = 30` is 0.668 and still climbing (0.26 → 0.42 → 0.68); the
> front simply has not uncovered the grain yet, and the recorded Jul-30 baseline never
> released either. Reading "no release by `t = 30`" as a failure — as rounds 1–3
> implicitly did — was wrong. Round 4 runs to **`t = 60`**, with a **matching `t = 60`
> control arm** so the comparison is like-for-like.

**Round-4 pass criteria** (revised accordingly): (1) `net = 0` — leak gone;
(2) **no depletion layer and no saturation** — `C_L(r)` around the grain flat and equal
to the far field, the round-3 and round-2 failure signatures respectively;
(3) `φ_liq` rising smoothly and at least as fast as the control, which no longer has a
grain eating its shell; (4) melt and enthalpy budgets unchanged; (5) tracer identity
closes. Jobs **20216530** (gate ✅ passed, `λ` bit-identical), **20216531** (control,
`t=60`), **20216532** (mask), **20216533** (mask + non-solid norm).

#### B.1.4-fix ROUND 4 RESULT — ✅ **the no-flux form works** (2026-08-29)

Jobs 20216530 (gate), 20216531 (control, `t=60`), 20216532 (mask), 20216533 (mask +
non-solid shell norm). All four to `t = 60`.

| | control | mask (Neumann) | mask + non-solid norm |
|---|---:|---:|---:|
| audit `net` | −4.405e-4 | **0.000000e+00** | **0.000000e+00** |
| `net` / `cum_melt` | **66 %** | **0** | **0** |
| `cum_melt` | 6.659e-4 | 6.566e-4 | 6.566e-4 |
| `max φ_liq` | 0.687 | 0.600 | 0.614 |

**Criterion 1 — leak gone.** `net = 0.000000e+00` exactly, as in round 3. Note the
control's leak has grown from 44 % of the melt budget at `t = 30` to **66 % at `t = 60`**:
it is not a bounded startup error, it accumulates for as long as the grain sits in
liquid, which is the whole production run.

**Criterion 2 — no depletion layer and no saturation.** `C_L` vs distance from the grain
surface at `t = 60`:

| `r − R` (cells) | −1:0 | 0:1 | 1:2 | 2:4 | 4:8 | 8:16 | 16:32 |
|---|---:|---:|---:|---:|---:|---:|---:|
| control | 0.657 | 0.693 | 0.688 | 0.708 | 0.696 | 0.680 | 0.641 |
| **round 4 (Neumann)** | 0.235 | **0.652** | **0.633** | **0.652** | **0.648** | **0.659** | **0.650** |
| round 3 (Dirichlet) | 0.040 | 0.033 | 0.100 | 0.255 | 0.440 | 0.536 | — |
| round 2 (rim restore) | — | — | — | — | — | — | ring saturated to 0.99 |

Flat from the grain surface outward and equal to the far field — the round-3 depletion
layer and the round-2 saturation are both absent. The `−1:0` band straddles the grain
surface, where the mask correctly reads ~0 and the **control counts liquid sitting inside
the rock**; that difference is the fix working, and it is also why the control's `φ_liq`
(0.687) is inflated above the mask's (0.600). The non-solid renormalisation recovers part
of it (0.600 → 0.614) by dividing out the `C_S` weighting of the partly-solid cells, as
designed — a real but modest effect, since `(1 − vf_cell)` already down-weights them.

**Criterion 4 — budgets.** `cum_melt` differs by 1.4 %, entirely accounted for by the
control melting extra ice into cells inside the grain.

> **⛔ `StageB_release2d` CANNOT test the release, and never could.** The control's
> `φ_liq` **plateaus and then reverses**: 0.6842 (`t=40`) → 0.6869 (`t=50`) → 0.6858
> (`t=55`) → 0.6836 (`t=60`). Its walls are no-flux, so the ocean's finite heat content
> is exhausted, the front stalls and slightly refreezes. No threshold near 0.9 is
> reachable on this deck at any run length. Rounds 1–3 read "no release" off this deck as
> a symptom; it is a property of the deck. **The release must be tested on
> `StageB_release_fn`**, which holds a Dirichlet warm bottom wall (`BC_BS/BC_CS = 1.0`)
> and is the deck that produced the recorded B.1.3-fn PASS (release at `t = 19.026013`,
> job 19581767).

**Release re-test submitted: jobs 20217237 (control) / 20217238 (mask + non-solid norm)**
on `ECCO_TESTS/StageB1_development/B14_relfn/`, both with `F_release = 0.9` to match the recorded baseline
rather than the 0.709 calibration, which was fitted to the old normalisation and is
exactly what is being re-tested. The control arm doubles as a harness check: it should
reproduce `t ≈ 19.03`.

#### B.1.4-fix ROUND 4 — RELEASE RE-TEST ✅ **PASS** (jobs 20217237 / 20217238, 2026-08-29)

Run on `StageB_release_fn` (warm Dirichlet bottom wall), the only deck that can actually
reach the threshold. Deck `ECCO_TESTS/StageB1_development/B14_relfn/`, `F_release = 0.9` in both arms to
match the recorded B.1.3-fn baseline.

| check | control | round-4 fix (no-flux + non-solid norm) | bar |
|---|---|---|---|
| **harness validity** | release at `t = 19.026013` — **identical to the recorded job 19581767** | — | must reproduce |
| audit `net` | −2.705e-4 (13.6 % of melt) | **0.000000e+00** | 0 |
| release fires | `t = 19.026013` | `t = 17.247013` | must fire |
| locked before release | `max\|U_y\| = 0.000e+00`, `\|ΔY\| = 0.000e+00` | `0.000e+00`, `0.000e+00` | exactly 0 |
| impulse-free handoff | 0.2593 = **0.52×** `a_free` | 0.2573 = **0.51×** `a_free` | < 2× |
| settles | `U_y = −1.004e-3`, final `y = 1.14394` | `U_y = −9.506e-3`→`−9.506e-4`, final `y = 1.14246` | falls |

The control arm reproducing the recorded release time **to every printed digit** is what
makes this comparison trustworthy: the deck, the binary lineage and the diagnostic are
all behaving as they did in July, so the only thing that moved is the physics under test.

**The lock, the ramp and the settling are untouched by the fix** — velocity exactly zero
for the whole locked phase, handoff peak 0.51× free-settling against 0.52× for the
control. B.1.3's mechanism survives the operator change intact.

**Why the release is 1.78 t.u. (9.4 %) earlier, and why that is the corrected value.**
Not the criterion: the recalibration below shows the renormalisation barely moves `φ`.
It is that the control's grain **eats liquid from its own surroundings**, which suppresses
melting near it — `cum_melt` is 1.990e-3 for the control against **2.193e-3 (+10 %)** with
the leak removed. The front therefore reaches any given exposure sooner once the grain
stops acting as a sink. (Job 20231742 runs the mask-only arm, old normalisation, to
confirm this attribution directly.)

##### Attribution of the 1.78 t.u. shift — ✅ **it is the physics, not the criterion** (job 20231742)

Third arm: CH no-flux mask ON, **old** shell normalisation, so the two changes separate.

| arm | binary | release | share of shift |
|---|---|---:|---:|
| control (leak present, old norm) | `stageB_sedctl` | `t = 19.026013` | — |
| **mask only** (leak removed, old norm) | `stageB_mask` | `t = 17.486013` | **1.540 t.u. = 87 %** |
| mask + non-solid norm | `stageB_mask2` | `t = 17.247013` | 0.239 t.u. = 13 % |

**87 % of the shift is the leak removal**, confirming the mechanism proposed above: the
control's grain acts as a liquid sink that suppresses melting in its own neighbourhood
(`cum_melt` 1.990e-3 → 2.192e-3, **+10 %**), so the front reaches any given exposure later.
Only 13 % is the criterion being renormalised — consistent with the recalibration table,
where the two normalisations agree to ~1 %. The corrected release time is therefore a
**physics correction, not a redefinition**, which is the answer needed for the write-up.

Both mask arms report `net = 0.000000e+00` and near-identical `cum_melt` (2.192493e-3 vs
2.192805e-3), as expected: the renormalisation touches only the release criterion, never
the solved fields.

##### `F_release` recalibration — ✅ **0.709 stands**

`φ_liq` against exposure (0 = front at the grain's underside, 1 = grain just uncovered),
measured on the round-4 fix run:

| exposure | 0.0 | 0.5 | **1.0** | 1.5 |
|---|---:|---:|---:|---:|
| `φ_liq`, **new** (non-solid) normalisation | 0.166 | 0.526 | **0.701** | 0.885 |
| `φ_liq`, old normalisation (B.1.3-fn recorded) | — | 0.565 | **0.709** | 0.867 |

The two agree to ~1 %, so **the B.1.3-fn calibration `F_release ≈ 0.709` carries over
unchanged** and still means "released the moment it is mechanically free". The
renormalisation's value is not a shifted calibration but the removal of the artificial
ceiling: `φ` is no longer capped by the grain's own `C_S` halo, so the
**`F_release ≥ 0.98` silent-hang trap of B.1.3-fn should be gone**. That last point is
*inferred, not directly verified* — this run released at 0.9 and the grain then fell away,
so `φ` stopped climbing at 0.945; a deck that never releases would be needed to confirm
the new ceiling reaches 1.0.

> **⚠ The `F_release = 0.9` in the B14_relfn decks is the BASELINE-MATCHING value, not the
> recommendation.** Production decks should use **0.709**. `0.9` fires ~0.67 diameters
> late, exactly as B.1.3-fn found.

#### B.1.4-fix — ✅ **RESOLVED. Round 4 is the answer.** (2026-08-29)

Four rounds, and the first three were wrong. Recorded in full above rather than
edited over, because the failure modes are the useful part.

| round | approach | outcome |
|---|---|---|
| 1 | restore into the `C_S`-threshold "rim" | ❌ the "rim" was a disc of 3.5 R (tanh tail); **melted the whole ice slab** |
| 2 | restore into a one-cell mesh dilation | ❌ ice survives, but the ring **saturates to `C_L≈0.99`** and releases spuriously |
| 3 | CH mask, **Dirichlet** `C_L = 0` in the grain | ✅ leak gone; ❌ **depletion layer 10 cells deep** — asserts the rock forces ice against itself |
| **4** | CH mask, **zero-FLUX** (Neumann) at the grain surface | ✅ **leak gone, profile flat, release mechanism intact, calibration unchanged** |

**Production flags for every ECCO deck with a resolved grain:**
`VOF_DIFFUSE_SEDIMENT_CH_MASK` + `VOF_DIFFUSE_SHELL_FRACTION_NONSOLID`.
`VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE` is **superseded — do not enable it**; with the mask
on there is nothing to restore, and on its own it produces the round-1/2 artifacts. Keep
`VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT` on in production: `net` is the standing proof that the
leak is zero, and it costs one reduction.

**What the leak actually was, now that it is measured rather than inferred.** Not a
bounded startup artifact and not 9.3 %: the signed `net` is **44 % of the melt budget at
`t = 30`, 66 % at `t = 60`** on the 2-D deck and 13.6 % by `t = 25` on the release deck,
growing for as long as the grain sits in liquid. The original 9.3 % figure came from the
meltwater-tracer identity, which conflates the leak with tracer transport error; the
`SED_CLIP_AUDIT` `net` column measures it directly.

**Non-invasiveness held throughout.** The Stage-A no-impact gate was re-run on every
round (jobs 20198502, 20213564, 20215987, 20216530) and returned
`λ = 0.2184495545645199` and enthalpy `4.975384548799866e-10` **bit-identical in all
16 digits every time**. None of these four rounds ever reached the validated Stage-A path.

**Three methodological corrections worth carrying forward:**
1. **Never key a locality argument on a tanh tail.** `C_S = 0.05` sits exactly on `r = R`,
   but the sub-cutoff set runs to `diffuse_solid_support_extra(Cn)` — 2.5 R beyond the
   surface here. Use the mesh (`diffuse_sediment_adjacent`), not a level-set threshold.
2. **A silent fallback hides a broken fix.** Round 1's residual was reported *after* the
   global band restore, so a restore that placed nothing locally still read
   `sed_resid ~ 1e-20`. Report what a local step actually placed, and drop rather than
   teleport what it could not.
3. **Check that the deck can express the failure you are testing for.** "Release never
   fires" was read as a symptom for three rounds; `StageB_release2d` has no sustained heat
   source, its `φ_liq` plateaus at 0.687 and *reverses*, and it can never reach 0.9 at any
   run length. Release must be tested on `StageB_release_fn`.

#### B.1.4-fix — ✅ **RESOLVED. Release re-test on `StageB_release_fn` passes** (2026-08-29)

Jobs 20217237 (control) / 20217238 (round-4 no-flux mask + non-solid shell norm), both
`F_release = 0.9`, both to `t = 25`.

| check | control | **round-4 fix** | bar |
|---|---:|---:|---|
| audit `net` | −2.705e-4 (**13.6 %** of melt) | **0.000000e+00** | 0 |
| release fires | `t = 19.026013` | `t = 17.247013` | must fire |
| **locked**: `max\|U_y\|` | 0.000e+00 | **0.000e+00** | exactly 0 |
| **locked**: drift `\|Δy\|` | 0.000e+00 | **0.000e+00** | exactly 0 |
| impulse-free: peak `\|dU_y/dt\|` | 0.2593 | **0.2573** (0.51× `a_free`) | < 2× |
| settles after release | −1.004e-3 | −9.506e-4 | falls |
| `cum_melt` | 1.990e-3 | 2.193e-3 (**+10 %**) | — |
| final `X₁` | 1.14394 | 1.14246 | same rest position |

**The control arm reproduces the recorded B.1.3-fn baseline exactly** — release at
`t = 19.026013`, peak `|dU_y/dt| = 0.2593`, identical to job 19581767 in every printed
digit. So the harness, deck and analyser are validated, and the differences in the fix
column are the physics change and nothing else.

**The earlier release is the correct direction, not a regression.** The control destroys
liquid in the ice/sediment overlap for the whole run, which keeps the grain's *own release
shell* artificially dry and therefore delays its release. Remove the leak and the shell
wets properly: the grain frees **9 % sooner** and **10 % more ice melts**. Both numbers
move the way conservation says they must, and the grain still settles to the same rest
position. Lock rigidity is still *exactly* zero for the whole locked phase, and the ramp
still delivers an impulse-free handoff at half the free-settling acceleration.

**Final form of the fix (all in `Boundary.h`, all `#undef` by default):**

| flag | what it does |
|---|---|
| `VOF_DIFFUSE_SEDIMENT_CH_MASK` | identity rows for sediment cells in the CH biharmonic **plus** `diffuse_sediment_noflux_correction`, which cancels the stencil face terms reaching into sediment. The grain is a hole in the CH domain with **no-flux** walls (90° contact angle — the rock is indifferent between water and ice). Requires `VOF_IBM`; `#error` under `AXISYM_RZ`. |
| `VOF_DIFFUSE_SHELL_FRACTION_NONSOLID` | `φ = ∫(1−vf)C_L / ∫(1−vf)(1−C_S)` — divides the grain's own `C_S` halo out of the release criterion |
| `VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT N` | log-only; `net = 0` is the standing proof the leak is gone |
| `VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE` | **superseded — do not use.** Round-2 rim restore; kept only so the round-1/2 result is reproducible |

**Stage-A no-impact gate passed on all four rounds** (20198502, 20213564, 20215987,
20216530): `λ = 0.2184495545645199` and enthalpy `4.975384548799866e-10`, bit-identical
every time. Nothing in this work reaches the Stage-A path.

> **✅ `F_release` recalibration DONE — 0.709 stands.** *(This block previously said
> recalibration was outstanding and that 0.945 was a new ceiling. The exposure sweep has
> now been run against the new definition and both statements were wrong; corrected here
> with the measurement.)*
>
> `φ_liq` vs exposure (0 = front at the grain's underside, 1 = grain just uncovered), from
> the round-4 fix run (job 20217238):
>
> | exposure | 0.0 | 0.5 | **1.0** | 1.5 |
> |---|---:|---:|---:|---:|
> | `φ_liq`, **new** (non-solid) normalisation | 0.166 | 0.526 | **0.701** | 0.885 |
> | `φ_liq`, old normalisation (B.1.3-fn recorded) | — | 0.565 | **0.709** | 0.867 |
>
> The two agree to ~1 %, so **`F_release ≈ 0.709` carries over unchanged** and still means
> "released the moment it is mechanically free". It is a **production value**, not a
> shakeout placeholder.
>
> **0.945 is not a ceiling.** It is simply where `φ` stopped climbing because the grain
> released at 0.9 and fell away from the ice. The renormalisation's actual effect is to
> *remove* the old `C_S`-halo cap, so the `F_release ≥ 0.98` silent-hang trap should be
> gone — **inferred, not directly verified**, since no run has yet let a fully-uncovered
> grain sit unreleased. Confirming the new ceiling reaches 1.0 is the only B.5 item left
> here.

#### B.1.4 Meltwater passive tracer (3rd scalar) — ✅ **IMPLEMENTED 2026-07-30**
- **`PARTIES/src/Eulerian/Conc.c`** — field 2 = meltwater tracer, `Ri = 0` (no buoyancy),
  injected at the interface at the local melt rate `V_Γ|∇F|`.
  Measures fresh-meltwater entrainment in the wake — the central ECCO quantity. **`NConc = 3`**.

`Conc_add_meltwater_RHS`, called from the RHS assembly for `iconc == 2` behind a new
`[phase_change] meltwater_tracer` flag (default **0 = off**, so decks without it are unchanged).

**It uses the *same discrete* `vof->melt_src` field that drives the phase change and pays the latent
heat** — deliberately, because that turns the tracer into a *conservation check* rather than a
decoration. With no sink and no-flux walls,

    d/dt ∫C_mw  =  ∫m  =  d/dt ∫F ,

so the tracer's total must track the ice volume lost, cell for cell; any drift between them is a
transport error in one of the operators. Refreezing (`m < 0`) correctly withdraws tracer, so the
balance holds through refreezing episodes too. `Ri = 0` for field 2 is required, not optional: the
actual freshening is already carried by the salinity field, so giving the tracer buoyancy would
double-count it. `[MPI]` pointwise, no stencil, rank-count independent.

#### B.1.5 Combined initial condition — ice layer at top + embedded sphere + stratified water
**No new VOF init type is needed — `init_type 28` already IS the ECCO configuration.**
`VoF_init_horizontal_ice_layer` sets `F = 0.5(1 - tanh((y - y_int)/(2√2 Cn)))`, i.e. **ice above,
water below**, with the interface height set by `vof_slab_x0·Ly` — exactly the ECCO geometry, and
already validated through the Favier melting-RB campaign. The earlier plan for a bespoke
`init_type = 29` carving a pocket for the grain was unnecessary: the grain is painted by the IBM as
`C_S`, and the ternary `C_L + C_S <= 1` clip makes it displace the ice automatically, so the ice
field needs no hole. (`init_type 29` therefore stays unused; note `28` is *not* free — it is the
Stage-A melting-RB layer, `Cart3d.c:1150`.)

**✅ IMPLEMENTED 2026-07-30** — what B.1.5 actually required was the two *scalar* inits, the
y-oriented analogues of the vertical-slab inits 30/31:

| `conc_init_type` | function | profile |
|---|---|---|
| **34** | `Conc_init_ice_layer_T` | `theta = theta_ice + (cbd0 - theta_ice)·Flayer(y)` — ice cold above, ocean `cbd0` below, blended with the *same* tanh as VOF init 28 so the scalar and phase fields are consistent at `t = 0` |
| **35** | `Conc_init_ice_layer_S_ylinear` | `s = [cbd5 + (cbd2 - cbd5)(y-ymin)/Ly]·Flayer(y)` — linearly stratified ocean, masked to zero inside the ice; `cbd5 > cbd2` gives the stable polar column |
| meltwater tracer | (default) | `Conc_init_zero` — no new code |

Both are computed analytically from the inputs rather than from `F`, so they do not depend on the
VOF/Conc initialization order. Unlike init 33 (Favier melting-RB) init 34 imposes **no** conductive
profile and **no** perturbation: the ECCO ambient is a stratified ocean, not an RB cell, and the
convection is driven by the melt itself.

- Particle: listed in **`p_mobile.inp`** with centroid inside the ice, and locked in place by the
  release mechanism until the melt front frees it (§B.1.3) — *not* `p_fixed.inp`.

#### B.1.6 Mixing diagnostics (postprocessing)
- **`PARTIES/src/IO/post_processing.c`** — horizontally-averaged `⟨w′S′⟩(y,t)`, `⟨w′T′⟩(y,t)`,
  effective diffusivity `κ_eff = −⟨w′c′⟩/∂_y⟨c⟩`; meltwater penetration depth; background-PE sorting
  for mixing efficiency `η = ΔBPE/(ΔBPE + ε·Δt)` (Winters 1995, on saved 3D fields). Cheap profile
  output every ~10 steps; full 3D fields via existing `Output_h5_data` + `ParticleOutput_h5`.
  **[MPI]** horizontal averages are `MPI_Allreduce` over the x–z pencil; the **BPE density sort must
  be a parallel/distributed sort** (or done offline on the HDF5 fields) — never an `MPI_Gather` of the
  whole field to rank 0, which would not scale at production grid sizes.

**✅ IMPLEMENTED OFFLINE 2026-07-30** — `ECCO_TESTS/StageB1_development/StageB_release2d/analyze_release.py`, taking the
roadmap's own "or done offline on the HDF5 fields" option. That choice **removes the distributed-sort
problem entirely** rather than solving it, which is the right trade: the sort is needed only at
output cadence, and doing it in numpy on saved fields costs nothing and cannot introduce a
rank-count dependence into the solver.

Implemented: `⟨w′c′⟩(y)`, `κ_eff = −⟨w′c′⟩/∂_y⟨c⟩`, meltwater penetration depth, background PE by
explicit sort, viscous dissipation `ε = (2/Re)S_ijS_ij`, and `η = ΔBPE/(ΔBPE + ε·Δt)`.

Two decisions worth recording:
- **Only liquid cells enter the BPE sort.** Ice is not fluid; including it would let the solid's
  buoyancy masquerade as available potential energy, and since the ice volume changes in time that
  would put a spurious trend straight into `ΔBPE`.
- **The buoyancy field used for sorting is stated, not assumed.** The shakeout is freshwater, so
  `b = θ`; a salty ECCO run must substitute the run's actual EOS. A wrong `b` makes BPE meaningless
  while still producing plausible-looking numbers, so this must never be left implicit.

Smoke-tested on the Stage-B freshwater run: `ΔBPE = +3.10e-3`, `ε = 1.84e-3`, `η = 0.144`.
Caveat carried in the output itself: `η` currently uses the final-time dissipation across the whole
span; a properly time-integrated `ε` needs the cheap profile output every ~10 steps that B.1.6
specifies, which is **not** yet implemented in-solver.

> **Checklist (Stage B code):** EOS re-fit (inputs only) ✅ (§B.1.1/B.1.1a — recommendation: use the
> **linear** buoyancy path, not `EOS_NONLINEAR`); `VOF_IBM` + `LAG_PARTICLE_RESOLVED` ✅ (B.0, both
> parts); interface-triggered **motion lock + release** ✅ (`PARTICLE_RELEASE` and its release list
> are *not* used, and no list transfer happens — §B.1.3); ECCO scalar inits 34/35 ✅ (no new VOF
> init needed — 28 already is the ECCO geometry); meltwater tracer (`NConc=3`) ✅; mixing
> diagnostics ✅ offline (in-solver profile output every ~10 steps still ☐). Collisions
> (`ACTM` + `LUBRICATION_NORMAL`) for the bottom wall reused as-is.
> Sediment-clip conservation ✅ **RESOLVED 2026-08-29** (§B.1.4-fix) — after three failed
> designs, by `VOF_DIFFUSE_SEDIMENT_CH_MASK` in its **no-flux** form plus
> `VOF_DIFFUSE_SHELL_FRACTION_NONSOLID`. Audit `net = 0`; release/lock/ramp all pass against the
> recorded B.1.3-fn baseline. `VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE` is superseded — do not use.
>
> **Remaining before the ECCO production campaign:** the §B.2 3-D shakeout (job 20218407, in
> flight) must pass; then `F_release` must be **recalibrated** against the non-solid shell
> normalisation (the 0.709 of B.1.3-fn was fitted to the old definition and the ceiling has moved
> 0.977 → 0.945); then the in-solver mixing-diagnostic profile output (B.1.6, still ☐) if the
> production campaign wants a time-integrated `ε`. The 2-D lock/ramp/release chain is **exercised
> and passing** on the resolved build (job 20217238); what remains unexercised is the same chain
> **in 3-D**, plus the wall collision.

### B.2 Debugging pass (reduced affordable 3D domain — code shake-out only)

**CLOSED 2026-09-08.** All checks pass, re-verified under the release-ramp fix in
20492853. The timestep comparison 20487690 FAILS and is **retired as an acceptance
item**: it targets the escape from the ice, which is an ill-conditioned
near-cancellation whose timing does not converge at any affordable dt. Accepted as
a documented constraint on B.3, not an open defect. See the closure report.
See [the closure audit](YANG_ECCO_B2_CLOSURE_REPORT.md) for the corrected physics,
new configuration, and current evidence. The runs below are retained history;
“in flight” and superseded pass claims in that history are not current status.

> **Run 1 verdict: physics unusable, infrastructure result valuable — do not discard the log.**
>
> It ran `parties.stageB_sed3d` built with the **round-1** sediment restore, which §B.1.4-fix
> ROUND 1 shows melts the ice slab. Its melt/release/settle behaviour means nothing. In
> particular the grain **never released** (still locked at `t = 125`), and the audit's healthy
> look (`sed_resid ~ 10⁻²⁰`, `ratio = 1.4 × 10⁻²`) was an artifact of the silent global fallback.
>
> **What it did establish, and what carries forward:**
> - The **3-D coupled path runs clean**: 12 h, 12 600+ steps, three scalars, IBM, HYPRE, 64 ranks,
>   **no NaN/Inf, no solver failure**, max velocity divergence steady at `3.5 × 10⁻¹⁵`. This was the
>   first time the full 3-D path had ever been run at all, and it did not fall over.
> - The binary's compiled configuration is confirmed correct from the run header: `x`,`z` periodic,
>   bottom/top no-slip, `TWOD_CARTESIAN` off.
> - **Walltime was ~2× short**: it reached `t = 125` of `time_max = 250` in 12 h. Run 2 needs
>   **24 h** (or a two-stage `Resume_h5` protocol).
>
> **Run 2 submitted 2026-08-29 as job 20218407**, from `t = 0` (the run-1 restart files carry a
> corrupted `C_L` field), 64 ranks, **24 h**, binary `parties.stageB_sed3d` rebuilt with the
> resolved B.1.4 fix (`..._CH_MASK` no-flux form + `..._SHELL_FRACTION_NONSOLID`).
> `F_release = 0.709` is kept for the shakeout, and per the recalibration above it is also the
> correct **production** value — "released the moment it is mechanically free".
>
> **Binary configuration** — `PERIODIC_NOSLIP_BOX` (x,z periodic; y no-slip: bottom wall, top
> anchors the ice) with `TWOD_CARTESIAN` **off**, plus the full Stage-A set, `VOF_IBM`,
> `LAG_PARTICLE_RESOLVED`, **`VOF_DIFFUSE_SEDIMENT_CH_MASK`** (§B.1.4-fix round 4 — the leak
> scales with sediment surface area, so 3-D is where it matters most; the superseded
> `..._CLIP_RESTORE` is NOT used), `VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT 200`,
> and **`LUBRICATION_NORMAL` on**. That last one is a deliberate change from the 2-D arms: B.2's
> check (4) is the bottom-wall collision, and exercising a collision path that production will not
> use would defeat the purpose. If it destabilizes, that is a bug B.2 exists to find.
>
> **Two deck corrections carried from the Jul-30 draft** (recorded in the deck's own header):
> `conc_init_type` is `{34, 35, 0}`, not the `{32, 33, 34}` written in the block below — the
> implemented ECCO inits are 34/35 and the tracer needs no init; and `Pe_T` is 100 with
> `time_max = 250`, not 700/20, because at `Pe_T = 700` the front needs `t ~ 1400` to clear the
> grain and the run would have ended with it still locked, looking like a release bug.
>
> **One deck change made now:** `F_release` `0.9 → 0.709`, the B.1.3-fn calibration. `0.9` fires
> 0.67 diameters late.

Purpose: exercise the **full 3D coupled path** (melt + release + settle + scalars + collisions)
cheaply, to find bugs — not to produce science. Reduced domain, reduced `Sc`, one particle, short
melt-out.

#### B.2 RUN 1 RESULT — ⚠️ **partial pass: 3 of 5 checks; the deck cannot test the other 2** (job 20218407, 2026-08-30)

Completed the full `time_max = 250` in 13 h 08 m on 64 ranks (faster than the 9.0 t.u./h
early estimate — `dt` grew to `max_dt`). Binary `parties.stageB_sed3d` with the round-4
no-flux CH mask.

| # | B.2 check | result |
|---|---|---|
| 1 | particle stationary while locked | ✅ **exact** — `X_y = 7.5` unchanged, `U = 0`, for all 25 000 steps |
| 2 | impulse-free release | ⛔ **NOT TESTED — never released** |
| 3 | meltwater tracer injected at the interface and conserved | ✅ **3.25e-3, and bounded** (3.25e-3 → 2.47e-3 → 3.25e-3 across `t = 50/125/250`) |
| 4 | bottom-wall collision | ⛔ **NOT TESTED — never released** |
| 5 | heat/salt/tracer budgets close | ✅ enthalpy drift **2.72e-4**, salt drift **0**, `C_S` drift **exactly 0** |
| — | B.1.4 leak in 3-D | ✅ `net = 0.000000e+00` — **the round-4 fix holds in 3-D** |

**What this establishes.** The 3-D coupled path is sound: no NaN, no solver failure, the
motion lock is exact over 25 000 steps, the sediment field is untouched, and the
meltwater tracer now closes to **3.25e-3 against the 9.3e-2 that opened B.1.4 — a 29×
improvement, in 3-D with a resolved grain, the hardest case.** Crucially the tracer error
is **bounded, not accumulating**, which is the qualitative difference from the old leak
(44 % → 66 % of the melt budget). It sits ~3× over the 1e-3 bar; the residual is at the
ice–water front, not the sediment (`net` is exactly zero), so it is a diffuse-interface
effect and not the B.1.4 defect returning.

**Why checks 2 and 4 could not run — the same defect as `StageB_release2d`.** Both `y`
walls were **no-flux for all three scalars** (`BC_A = 1, BC_B = 0, BC_C = 0`), making the
ocean a **finite heat reservoir**. Sensible heat available is `∫θ = 0.696`; melting all
the ice costs `St⁻¹ × 0.304 = 4 × 0.304 = 1.217`. **At most ~57 % of the ice can ever
melt.** Measured: `mean(θ)` decayed 0.696 → 0.405 while liquid fraction rose 0.6957 →
0.7685, and `φ_liq` climbed to 0.565 (threshold 0.709) with per-25-t.u. increments
decaying monotonically 0.058 → 0.021. Reaching the threshold would need several hundred
more time units against a shrinking driving temperature — not worth buying.

> **This is the second time this bit us, and the lesson from round 4 was not applied
> widely enough.** The round-4 note said "check that the deck can express the failure you
> are testing for" — but only `StageB_release2d` was checked. `StageB_release3d` carries
> the identical no-flux thermal BC and was submitted anyway.

**Fix and rerun.** `StageB_release3d/parties.inp` now takes a **Dirichlet `θ = 1` on the
south (bottom) wall**, matching `StageB_release_fn`, the deck that demonstrably releases;
top wall and both scalar fields stay no-flux. Previous deck kept as
`parties.inp.noflux_bak`. **Rerun submitted as job 20240184.**

> **⚠ §B.3 PRODUCTION HAS THE SAME DEFECT.** The B.3 input block below specifies
> "no-flux T,S,tracer on bottom+top walls (`BC_A[NS]=0, BC_B[NS]=1, BC_C[NS]=0`)". A
> production melt-out-and-settle run under that BC will stall the same way — the box can
> only melt the ice its own initial heat content pays for. Either give it a warm bottom
> wall, or adopt the two-stage `Resume_h5` protocol the sweep paragraph already
> contemplates (stage A melt-out, stage B settling restart). **Decide before sizing the
> campaign.**

#### B.2 RUN 2 — warm bottom wall: ✅ **melting is now sustained**, but the deck is still too short (job 20240184, 2026-08-30)

Full `time_max = 250` in 15 h 02 m. Same three checks pass as run 1 (lock exact,
`net = 0.000000e+00`, budgets closed, no NaN). The BC fix did what it was meant to:

| | run 1 (no-flux) | **run 2 (warm bottom wall)** |
|---|---|---|
| `mean(θ)` at `t = 0 / 125 / 250` | 0.696 / 0.518 / **0.405** (collapsing) | 0.696 / 0.578 / **0.546** (**stabilising**) |
| `φ_liq` increment per 25 t.u. | 0.058 → **0.021**, monotone decay | 0.058 → **0.038**, plateaued |
| `φ_liq` at `t = 250` | 0.565 | **0.644** |
| release | never, and never would | not yet — but on a steady approach |

**The diagnosis is confirmed and the fix is correct.** Under no-flux the reservoir
collapsed and the front was asymptoting *below* the threshold; with the Dirichlet wall
`mean(θ)` levels off near 0.546 and `φ_liq` rises at a **constant** ~0.038 per 25 t.u.
The deck was simply not long enough: 0.6436 → 0.709 at that rate needs **~43 more time
units, i.e. release near `t ≈ 293`**, with more beyond that for the settle and the
bottom-wall collision that checks (2) and (4) require.

**Resumed rather than rerun** — job **20253869**, restarting from the `t = 250` checkpoint
with `time_max = 450`. That buys the release plus ~157 t.u. of settling for ~12 h, against
15 h to redo what is already done. `jobscript_resume.sh` now takes `PREV` and `TIME_MAX`
as environment overrides so this is repeatable.

##### ⚠ Latent bug found in the restart path — `Resume.h5` alone is NOT sufficient (job 20253869, 2026-08-30)

The first resume attempt **segfaulted on startup**:

```
!!! Open dataset "/u" from @)<garbage> failed !!!
Caught signal 11 (Segmentation fault: address not mapped to object)
```

**Cause.** `src/IO/Resume.c` has **two** read paths and they open **different files**:

| line | opens | holds |
|---|---|---|
| `Resume.c:67` | `Resume.h5` | restart metadata / state |
| `Resume.c:382` | `Data_<abs(params->noutput)>.h5` | the **Eulerian fields** (`/u`, …) |

Staging only `Resume.h5` — the obvious thing to do, and what the roadmap's own restart
discussion implies — leaves the second `H5Fopen` returning a negative `file_id`. **That
return is not checked**, so the garbage handle is passed to `H5Dopen`, which prints the
uninitialised filename and dies with SIGSEGV rather than a clean "file not found".

**Fix in the jobscript**: carry the last two `Data_*.h5` snapshots alongside `Resume.h5`
and the `Particle_*.h5` files. `StageB_release3d/jobscript_resume.sh` now does this and
takes `PREV` / `TIME_MAX` as environment overrides.

> **This matters well beyond the shakeout.** §B.3 explicitly plans a *"two-stage protocol
> (`Resume_h5`) if the melt and settling timescales separate too far: stage A melt-out,
> stage B settling restart"* — and given B.2 has now needed a restart twice, that protocol
> is likely to be used, not optional. Anyone driving it must stage the `Data` snapshot too.
> The unchecked `H5Fopen` at `Resume.c:385` is worth fixing properly (a clean abort with
> the attempted filename) so the next person loses a queue slot instead of an afternoon.

#### B.2 RUN 3 (resume to `t = 450`) — ✅ **RELEASE FIRES IN 3-D**, and it exposed a penalization bug (job 20253930, 2026-08-31)

9 h 14 m from the `t = 250` checkpoint. **The grain released at `t = 311.378`**, `φ_liq`
crossing the calibrated `F_release = 0.709` threshold exactly.

| # | B.2 check | result |
|---|---|---|
| 1 | stationary while locked | ✅ **exact** — `max\|U\| = 0.000e+00`, `\|ΔY\| = 0.000e+00` over 31 138 steps |
| 2 | impulse-free release | ✅ peak `\|dU_y/dt\| = 0.1457` = **0.24×** `a_free` (bar < 2×) |
| 3 | tracer conserved | ✅ (run 1 measurement stands) |
| 4 | settles + wall collision | ❌ **grain creeps, never reaches the wall — see below** |
| 5 | budgets | ✅ `net = 0.000000e+00`, no NaN |

**The lock and the release are validated in 3-D.** Combined with B.1.3-fn in 2-D, the
mechanism is now exercised in both.

##### ⛔ BUG — the Brinkman ice mask counts the resolved grain as ice

Post-release the grain falls at a **constant** `U_y ≈ −5.7e-3` from the first step, never
accelerating (`t = 312 / 320 / 350 / 400 / 450` → `−5.54 / −5.55 / −5.63 / −5.74 / −5.71
e-3`). At `Ga = √((ρ_s/ρ_f−1)g d³)/ν = 122` the terminal velocity should be **O(1)** — the
grain is ~250× too slow, and it had only fallen 0.84 diameters by `t = 450`.

**Cause.** The `ICE_PENALIZATION` face mask is built as `φ_s = 1 − C_L`:

| file | line |
|---|---|
| `msolve_cg.c` | `123` (`vel_jacobi_diag`), `343` (`matVec`) |
| `msolve_direct.c` | `271` (`Velocity_solve_explicit`) |

Inside a resolved sediment grain `C_L = 0`, so `φ_s = 1` and the grain receives the **full
Darcy damping at `darcy_tau = 1e-4`** — on top of the IBM forcing that already represents
it. The rock is being treated as ice and pinned. This is precisely the round-1 "grain
creeps while locked" signature, but now acting *after* release, where it silently
suppresses settling instead of producing an obvious lock violation.

**Fix.** The ice fraction is `1 − C_L − C_S`; the resolved solid is not ice. All three
sites now subtract the face-averaged `C_S`, with `fsol == NULL` reproducing the previous
mask exactly. **Non-invasive by construction**: `C_S ≡ 0` in a Stage-A build, so the
expression is unchanged there.

> **Note this was invisible to every 2-D test.** In `StageB_release2d` the front stalls
> and the grain never releases; in `B14_relfn` the grain does release, but its recorded
> settling (`U_y ≈ −1.0e-3`, matching the July baseline **exactly**) was *also* damped by
> this same bug — the control and fix arms agreed with each other and with the recorded
> value, so nothing looked wrong. **A reproducible number is not a correct one.** It took
> a 3-D run with a known analytic expectation (`Ga = 122` ⇒ `U_t = O(1)`) to expose it.
> The B.1.3-fn "settles after release ✅" verdict should be read as "moves downward",
> not as a validated settling velocity.

**Verification submitted:** job **20268028** (Stage-A no-impact gate — must stay
bit-identical) and job **20268029** (resume from the `t = 450` state, where the grain is
already free at `y = 6.66`, to confirm it now accelerates to `O(1)` and reaches the wall).
A clean full rerun from `t = 0` is held until those pass.

##### B.2 check (4) — ⛔ **settling is broken, and the penalization mask was only part of it**

**Verification of the mask fix (jobs 20268028 / 20268029, 2026-08-31):**

| | result |
|---|---|
| Stage-A no-impact gate | ✅ `λ = 0.2184495545645199`, enthalpy `4.975384548799866e-10` — **bit-identical, 5th consecutive round** |
| settling after the fix | `U_y`: −5.71e-3 → **−7.36e-3** (+29 %) — real, but still **~190× too slow** |

So `φ_s = 1 − C_L − C_S` was a genuine bug and is genuinely fixed, but it is **not** the
main cause. Two measurements at `t = 480` locate the rest:

1. **The grain is no longer in ice.** Ice fraction in its shell is 0.014–0.025 and `C_L`
   is 0.955–0.988 — essentially pure water. Darcy damping cannot be holding it.
2. **The IBM is not enforcing rigid-body motion.** Mean `v_y` *inside* the grain
   (`r < 0.6R`) is **−2.28e-3** against a particle velocity of **−7.37e-3** — the fluid
   the particle is supposed to be carrying is moving at **31 %** of it. Meanwhile
   `max|v| = 0.79` elsewhere, so the solver is not stalled; only the particle coupling is
   weak.

That points at the IBM coupling — `N_forcing_loops = 1` at `d/Δx = 16` — rather than at
anything in the ECCO machinery.

> **⚠ NOTHING IN THIS PROJECT HAS EVER VALIDATED A SETTLING VELOCITY.** B.1.3-fn recorded
> "settles after release ✅", but that check only tests `U_y < 0`. Its measured
> `U_y ≈ −1.0e-3` reproduced the July baseline **exactly**, and the B.1.4 control and fix
> arms agreed with it and with each other — so three runs sharing one defect looked like
> confirmation. **Agreement between runs is not validation against physics.** The 3-D case
> only exposed it because it has an independent analytic expectation
> (`Ga = √((ρ_s/ρ_f−1)g d³)/ν = 122` ⇒ `U_t = O(1)`).

**New gate: `ECCO_TESTS/StageB2_closure/Settle_validation/`** — one sphere, quiescent liquid, **no ice, no
phase change, no release mechanism** (`stefan = 0`, `vof_slab_x0 = 0.999`,
`F_release = −1`), settling from rest. It isolates the IBM from every ECCO term and
measures `U_t` against the `Ga = 122` expectation. Two arms differing only in
`N_forcing_loops` (**1** vs **5**, jobs **20269208** / **20269209**) to test the leading
hypothesis directly.

**This gate should have existed before B.2.** Settling and the wall collision are half of
what §B.2 is for, and the resolved-particle machinery was inherited as "reuse as-is"
(§B.1.2) without a quantitative check. §B.3's whole science output is a settling
trajectory and its wake — that number cannot rest on an unvalidated terminal velocity.

##### SETTLING GATE RESULT — ⛔ **settling is broken in the simplest possible case** (jobs 20269208 / 20269209, 2026-08-31)

One sphere, quiescent liquid, **no ice, no phase change, no release mechanism**, from rest.
`Ga = 122` ⇒ `U_t = O(1)`. Free fall alone gives `a_free·t = 0.6t`.

| `t` | 1 | 5 | 10 | 20 | 40 | free fall at `t=40` |
|---|---:|---:|---:|---:|---:|---:|
| `U_y`, `N_forcing_loops = 1` | −7.27e-3 | −6.82e-3 | −6.59e-3 | −6.31e-3 | **−6.13e-3** | −24.0 |
| `U_y`, `N_forcing_loops = 5` | −2.49e-3 | −2.35e-3 | −2.25e-3 | −2.25e-3 | **−2.15e-3** | −24.0 |

**200–500× too slow, and the velocity is flat from `t ≈ 1` — the particle never
accelerates at all.** In 40 time units it falls 0.27 d (`nfl=1`) or 0.10 d (`nfl=5`).

**More IBM forcing makes it slower**, which rules out under-converged forcing and points
the other way: the particle is being *held*.

**Cause — a phantom ice field, amplified by an over-stiff `darcy_tau`.** Measured in this
**ice-free** deck at `t = 40`, where `φ_s = 1 − C_L − C_S` should be ≈ 0 everywhere:

| region | `φ_s` mean | Darcy rate `2φ_s/τ` |
|---|---:|---:|
| `r < 0.5R` | 0.171 | 3,410 |
| **`0.5R – R`** | **0.773** | **15,467** |
| `R … R+1` cell | 0.023 | 465 |
| `R+1 … R+3` | 0.040 | 802 |
| far field `r > 2` | 0.004 | 88 |

against an inertial scale `ρ/(α_k Δt) ≈ 227`. The Brinkman term beats inertia by **68×**
in the grain's diffuse rim.

Two compounding causes, both of which must be fixed:

1. **The ternary deficit is read as ice.** `C_L + C_S < 1` in the grain's diffuse solid
   rim, and the shortfall is charged to the ice phase. **The round-4 CH mask makes this
   worse by construction**: pinning `C_L = 0` wherever `C_S ≥ 0.05` (right for the CH
   problem) leaves `φ_s = 1 − C_S`, which is ~0.95 at the cutoff. *This is an interaction
   between the B.1.4 fix and `ICE_PENALIZATION` that was not anticipated when the mask was
   designed.* The `1 − C_L − C_S` mask correction of B.2 run 3 was necessary but not
   sufficient — the Brinkman term must be switched **off inside the resolved solid's
   support entirely**, since the IBM already governs there.
2. **`darcy_tau = 1e-4` is 3,906× stiffer than §A.1.5's own sizing rule.** That rule is
   `darcy_tau ≈ Re·Δx²`, chosen so the Brinkman layer `√(ν·τ)` is ~1 cell. For the 3-D
   grid (`Δx = 0.0625`, `Re = 100`) that is **0.391**; the deck's `1e-4` gives a layer of
   **0.016 cells**. The value was inherited from the 2-D Yang decks and never re-sized for
   this grid. Over-stiffness turns any spurious `φ_s` into a clamp.

**Third arm submitted: job 20273250** — identical settling deck with
`darcy_tau = 0.391` per the roadmap rule, to separate cause (2) from cause (1).


##### ROOT CAUSE — ✅ **`darcy_tau` was 3,906× too stiff. Settling now works.** (job 20273250, 2026-08-31)

One input, changed to the value §A.1.5 already prescribes:

| deck | terminal `U_y` |
|---|---:|
| `darcy_tau = 1e-4` (as written in every ECCO deck) | **−6.1e-3** |
| **`darcy_tau = 0.391`** (`≈ Re·Δx²`, the A.1.5 rule) | **−0.56** |

**A 90× change from a single input.** The particle now accelerates from rest and reaches a
clean terminal velocity by `t ≈ 5`, instead of sitting at a constant creep from `t ≈ 1`.

§A.1.5 states the rule and the reason: *"A useful first `darcy_tau` estimate is
`darcy_tau ≈ Re·Δx²`, which gives a Brinkman damping layer `√((1/Re)·darcy_tau)` of order
one grid cell."* The decks use `1e-4`, giving a layer of **0.016 cells**.

| deck | `Δx` | `Re` | rule value | deck `1e-4` is |
|---|---:|---:|---:|---:|
| 3-D shakeout `d/Δx = 16` | 0.0625 | 100 | **0.391** | 3,906× too stiff |
| **B.3 production `d/Δx = 24`** | 0.0417 | 100 | **0.174** | **1,736× too stiff** |
| B.3 hero `d/Δx = 48` | 0.0208 | 100 | 0.043 | 434× too stiff |
| 2-D `release2d` / `relfn` | 0.0039 | 500 | 0.0076 | 76× too stiff |
| Yang Stage-A `N = 288` | 0.0035 | 1000 | 0.012 | 121× too stiff |
| Yang Stage-A `N = 1440` | 0.00069 | 1000 | 0.00048 | 5× too stiff |

> **✅ Stage A is NOT affected and does not need revisiting.** With no resolved particle,
> `φ_s = 1 − C_L` is the true ice fraction and an over-stiff `τ` simply makes the ice more
> rigid — conservative, which is what that term is for. The `τ`-sweep null recorded in
> §G.8 ("τ = 2.5e-4 closed 0.8 % of the gap") is consistent with this: in the Stage-A
> regime the answer is insensitive to `τ`. The value only becomes destructive once a
> **resolved particle** is present, because then any spurious `φ_s` near the grain is
> multiplied by `2/τ`.

##### Second fix — Brinkman is now switched OFF inside the resolved solid

Independently of `τ`, `φ_s = 1 − C_L − C_S` was **0.77 over `0.5R–R` in an ice-free deck**,
because the ternary deficit (`C_L + C_S < 1` in the grain's diffuse rim) is charged to the
ice phase. **The round-4 CH mask makes this worse by construction** — pinning `C_L = 0`
wherever `C_S ≥ 0.05` leaves `φ_s = 1 − C_S`, ≈0.95 at the cutoff. An unanticipated
interaction between the B.1.4 fix and `ICE_PENALIZATION`.

All three mask sites (`msolve_cg.c` `vel_jacobi_diag` + `matVec`, `msolve_direct.c`
`Velocity_solve_explicit`) now set `φ_face = 0` where `C_S ≥ DIFFUSE_SOLID_MASS_CUTOFF`:
the IBM governs inside the resolved solid, so the Brinkman term must not act there.
`DIFFUSE_SOLID_MASS_CUTOFF` was promoted from a private `#define` in `VOF_DIFFUSE.c` to
`Eulerian/Include/VOF_DIFFUSE.h` so the two definitions cannot drift.

**Verification in flight:** job **20273751** (Stage-A no-impact gate), **20273752**
(settling, fix + `τ = 1e-4` — isolates the fix alone), **20273753** (settling, fix +
`τ = 0.391` — the production candidate).

##### ☐ STILL OPEN — the terminal velocity is right in magnitude but not yet validated

`U_t = 0.56` gives `Re_p = 56`; Schiller–Naumann for an unbounded sphere at that `Re_p`
predicts `Cd = 1.45` ⇒ **`U_t = 1.17`**. The measurement is **0.48× the correlation**, and
blockage is only 4.9 % of the cross-section — nowhere near enough to explain a factor of 2.

At `d/Δx = 16` a coarse IBM plausibly over-predicts drag, and production specifies
`d/Δx ≥ 24`, but **this factor of 2 is not yet accounted for**. A settling number that is
merely the right order of magnitude is not good enough for §B.3, whose entire science
output is a settling trajectory and its wake.

**Required before B.3: a resolution sweep** `d/Δx ∈ {16, 24, 32}` on the standalone
settling gate, against the drag correlation, with the confinement correction stated. These
runs are ~1 SU each — the cheapest remaining item on the critical path, and the last one
between here and a defensible production campaign.


##### SETTLING GATE — ✅ **RESOLVED. Both fixes were needed; the wall collision works.** (jobs 20273250 / 20273752 / 20273753)

Standalone gate: one sphere, quiescent liquid, no ice, no phase change, no release.
`Ga = 122`.

| arm | `darcy_tau` | Brinkman-in-solid fix | terminal `U_y` (mid-fall) | reaches the wall |
|---|---|---|---:|---|
| baseline | 1e-4 | no | **−0.0061** | no (0.27 d in 40 t.u.) |
| phantom-ice fix only | 1e-4 | yes | **−0.0112** | no |
| `τ` re-sized only | 0.391 | no | **−0.5628** | ✅ `t = 14.1` |
| **both fixes** | **0.391** | **yes** | **−1.0005** | ✅ **`t = 8.6`** |

**`darcy_tau` is the dominant cause (92× on its own); the Brinkman-in-solid fix nearly
doubles what remains (0.563 → 1.000). Both are required.**

**✅ B.2 check (4) — wall collision — is demonstrated here.** Both `τ`-corrected arms land
at `y = 0.4996` (sphere of `R = 0.5` resting on the wall) and come to rest at
`U_y ≈ −1.5e-5`, with no blow-up and no NaN. The collision path (`ACTM` +
`LUBRICATION_NORMAL`) works.

**Against the drag correlation.** `U_t = 1.000` ⇒ `Re_p = 100` ⇒ Schiller–Naumann
`Cd = 1.09` ⇒ self-consistent `U_t = 1.35`. The measurement is **0.74×** the unbounded
correlation, up from 0.48× before the second fix. The residual 26 % is attributable to
confinement (4d × 4d periodic cross-section with a return flow, plus the bottom wall) and
to `d/Δx = 16`; production specifies `d/Δx ≥ 24`. **This is now a quantification task, not
a bug** — the resolution sweep `d/Δx ∈ {16, 24, 32}` remains open, but the pathology is
gone and the magnitude is physical.

##### ✅ Do the Stage-B changes affect Yang / melting-RB? **NO — verified by measurement** (job 20274565)

This mattered because **the Le=100 no-impact gate has max speed exactly 0** and therefore
cannot exercise the velocity solver or `ICE_PENALIZATION` at all — while
`ICE_PENALIZATION` **is** defined in the Stage-A/Yang/RB flag set, and the 2026-08-31 work
edited `msolve_cg.c`, `msolve_direct.c` and `Velocity.c`. Six bit-identical gate passes did
not cover this. So a **convective** check was run: Yang freshwater `N = 288` to `t = 10`
against the recorded reference (`YangFresh_07272026/N288`), ice volume at matching output
times:

| `t` | 2 | 4 | 6 | 8 | 10 |
|---|---:|---:|---:|---:|---:|
| rel. diff | 3.96e-12 | 4.03e-12 | **4.06e-12** | 2.30e-09 | 1.36e-08 |

Against this project's own accepted bar — the Stage-A binary lineage test recorded
**5.8e-5 at `t = 6`** against a **7.2e-3** discriminating size — this is **seven orders of
magnitude tighter**. The growth to 1e-8 by `t = 10` is the documented Lyapunov
amplification (0.652 per time unit) of last-bit differences, exactly the B.0 part-2
signature.

**Conclusion: Yang and melting-RB do NOT need recomputing.** The mechanism is that
`C_S ≡ 0` without particles, so `cs_face = 0.0` exactly and `x − 0.0 == x` in IEEE, and the
`cs_face ≥ cutoff` branch never fires — now confirmed in a flowing case rather than
asserted. Melting-RB was not run directly (it compiles `PERIODIC_NOSLIP_BOX` instead of
`PERIODIC_Z_NOSLIP_BOX`), but the changed code and the `C_S ≡ 0` argument are identical and
Yang is the stricter test, being chaotic.

> **Separately: `darcy_tau = 1e-4` is also 121× stiffer than the A.1.5 rule for the Yang
> `N = 288` deck.** This does **not** invalidate any Yang result — with no particle,
> over-stiffness only makes the ice more rigid, which is the conservative direction and is
> consistent with the §G.8 `τ`-sweep null. Worth re-sizing if Yang is ever re-run for other
> reasons, but it is not a defect in the recorded results.


##### Settling resolution sweep + clean B.2 rerun — submitted 2026-09-01

**Sweep** (`ECCO_TESTS/StageB2_closure/Settle_validation/res{24,32}`, jobs **20397138** / **20397139**).
`Cn = 0.75Δx` and `darcy_tau = Re·Δx²` **both** scale with the mesh, so all three are
re-derived per arm — changing `NXM` alone would alter the Brinkman layer thickness and the
CH band and confound the convergence being measured:

| `d/Δx` | grid | cells | `Cn` | `darcy_tau` | ranks | `U_t` |
|---:|---|---:|---:|---:|---:|---:|
| 16 | 64×160×64 | 0.66 M | 0.046875 | 0.391 | 64 | **−1.0005** (measured) |
| 24 | 96×240×96 | 2.21 M | 0.031250 | 0.17361 | 64 | — |
| 32 | 128×320×128 | 5.24 M | 0.023438 | 0.09766 | 128 | — |

`time_max = 15` — the grain reaches the wall at `t ≈ 8.6` at `d/Δx = 16` and sooner at
finer resolution, so 40 was wasteful. Target: Schiller–Naumann `U_t = 1.35` unbounded at
`Re_p ≈ 100`; `d/Δx = 16` currently gives **0.74×**, and the sweep separates mesh
convergence from the 4d×4d periodic confinement.

**Clean B.2 rerun** (job **20397146**, `time_max = 400`, 36 h, 64 ranks). First run from
`t = 0` carrying every fix at once:

| fix | why runs 1–3 could not have passed check (4) |
|---|---|
| `darcy_tau` 1e-4 → **0.391** | 1e-4 clamped the freed grain to `U_y = 6e-3` against `U_t = 1.0` |
| Brinkman off inside the resolved solid | phantom ice field, `φ_s = 0.77` in the grain's rim |
| warm Dirichlet bottom wall | no-flux made the ocean a finite reservoir; the grain never freed |
| `time_max` 250 → **400** | release lands near `t = 311`, then ~7 d of settling to the wall |

`F_release = 0.709`, `Cn = 0.046875` (`= 0.75Δx`), `conc_init_type = {34,35,0}`,
`Pe = {100,700,100}` — all verified in the deck before submission. The standalone gate
(20273753) already shows `U_t = −1.0005` and a clean wall landing at `y = 0.4996`; this run
confirms it in the **coupled** case, where the grain must melt free first.

Expect it to run **faster per time unit than run 2's 15 h/250 t.u.**: the Brinkman diagonal
is 3,906× less stiff, so the velocity Helmholtz solve should converge in fewer iterations.


##### SETTLING RESOLUTION SWEEP — ✅ **CONVERGED, and the deficit is confinement** (jobs 20397138 / 20397139)

| `d/Δx` | grid | cells | `U_t` (mid-fall) |
|---:|---|---:|---:|
| 16 | 64×160×64 | 0.66 M | −1.0005 |
| 24 | 96×240×96 | 2.21 M | **−0.7724** |
| 32 | 128×320×128 | 5.24 M | **−0.7610** |

**Converged**: 24 → 32 moves `U_t` by **1.5 %**. Note `d/Δx = 16` *over*-predicts and the
sequence converges downward, so the earlier 1.0005 was an under-resolution artifact.

**The converged value is correct once confinement is accounted for:**

| quantity | value |
|---|---:|
| converged `U_t` | 0.761 (`Re_p = 76`) |
| Schiller–Naumann unbounded at that `Re_p` | 1.268 (`Cd = 1.243`) |
| measured / unbounded | **0.600** |
| Faxén wall factor at `d/D = 0.25` (Stokes limit ⇒ *maximum* retardation) | 0.506 |
| measured / Faxén | **1.19×** |

19 % less retarded than the Stokes limit — the correct direction at `Re_p = 76`, where
inertia weakens wall effects. **The 40 % deficit is the 4d × 4d periodic confinement, not
a solver error.** `d/Δx ≥ 24` is confirmed adequate, as §B.3 already specifies.

> The `d/Δx = 32` arm shows TIMEOUT at its 8 h limit but is **complete for its purpose**:
> it reached `t = 12.8` and the grain had already landed (`y = 0.4996`).

##### B.2 CLEAN RERUN — ✅ **all five checks pass**, ⛔ **but on a non-rigid ice slab** (job 20397146)

`t = 0 → 400` in 17 h 36 m, 64 ranks.

| # | check | result |
|---|---|---|
| 1 | stationary while locked | ✅ **exact** — `max\|U_y\| = 0.000e+00`, `\|ΔY\| = 0.000e+00` |
| 2 | impulse-free release | ✅ `0.4073` = **0.68×** `a_free` (bar < 2×) |
| 3 | tracer | ✅ `net = 0.000000e+00`, no NaN |
| 4 | **settles + wall collision** | ✅ **`U_t = −0.5827`, lands at `y = 0.4996`, comes to rest at `U_y = −1.1e-6`** |
| 5 | budgets | ✅ |

**Checks (2) and (4) are satisfied in the coupled case for the first time.** The grain
melts free, releases impulse-free, settles ~7 d, and lands on the wall without blow-up.

**But `darcy_tau = 0.391` is too soft for the ice.** Measured at `t = 30` in bulk ice
(`1 − C_L − C_S > 0.95`):

| | bulk-ice `max\|u\|` | `mean\|u\|` | as % of peak domain flow |
|---|---:|---:|---:|
| `τ = 1e-4` (run 2) | 1.6e-5 | 1.4e-6 | **0.00 %** |
| **`τ = 0.391`** | **0.277** | **0.096** | **14.1 %** |

**The ice slab is flowing.** That is what moved the release from `t = 311` (run 3) to
`t = 55.4` here — a mobile "ice" carries heat to the front, melting it 5.6× faster. **The
melt and release *timing* in this run are not trustworthy**; the particle mechanics
(lock, ramp, settling, collision) are, since they do not depend on ice rigidity.

**Root of the mistake:** §A.1.5 prescribes *"run a sensitivity sweep and choose the largest
damping time that makes ice velocity leakage irrelevant"*. That sweep was listed as open
item 3 and then **skipped** — `0.391` was taken straight from the 1-cell-layer rule, which
is a *starting estimate*, not the answer. The two requirements pull opposite ways:

| `τ` | bulk-ice rate `2/τ` | near-particle rate `2×0.02/τ` | vs inertia ≈ 227 |
|---:|---:|---:|---|
| 1e-4 | 20,000 | 400 | ice rigid ✅, **particle clamped** ❌ |
| **1e-3** | **2,000** | **40** | rigid (9×) ✅, mild (18 %) ✅ |
| 3e-3 | 667 | 13 | marginal (3×), free ✅ |
| 0.391 | 5 | 0.1 | **ice flows** ❌, particle free ✅ |

**`τ` sweep submitted on BOTH diagnostics** — particle (standalone settling gate, no ice)
and ice (short 3-D melt to `t = 35`, bulk-ice `max\|u\|`) — at `τ = 1e-3` (jobs 20421162 /
20421163) and `τ = 3e-3` (20421164 / 20421165). Pass = ice leakage ≲ 1 % of peak flow
**and** `U_t` within a few % of the 0.761 converged value.

> **B.2 is therefore NOT yet closed.** Four of five checks are resolution-independent and
> pass; the melt/release timing needs re-running at the `τ` this sweep selects. That is one
> more ~1,600 SU run, not a redesign.


##### `darcy_tau` SWEEP — ⛔ **there is NO viable τ. The two requirements are 2 orders apart.** (jobs 20421162–20421165, 2026-09-02)

Both diagnostics, at each τ: **particle** = standalone settling gate (no ice, `U_t`
target 0.761 converged / ~1.00 at `d/Δx = 16`); **ice** = short 3-D melt, bulk-ice
`max|u|` as a fraction of peak flow.

| `τ` | ice leakage | particle `U_t` | verdict |
|---:|---:|---:|---|
| 1e-4 | 0.00 % ✅ | −0.011 ❌ | ice rigid, particle clamped |
| 1e-3 | 0.02 % ✅ | −0.028 ❌ | ice rigid, particle clamped |
| 3e-3 | 0.19 % ✅ | −0.043 ❌ | ice rigid, particle clamped |
| 0.391 | **14.06 %** ❌ | −1.000 ✅ | **ice flows**, particle free |

Ice rigidity needs `τ ≲ 3e-3`; the particle needs `τ ~ 0.4`. **No overlap.**

**Why — a feedback, not a tuning problem.** `φ_s` measured in the shell around the grain
is **0.15–0.30, not the ~0.02 assumed when the window was proposed**, and it is *larger*
at weaker damping (0.161 at `τ=1e-3`, 0.302 at `τ=3e-3`). Brinkman damps the flow near the
grain ⇒ CH cannot relax `C_L` up to `1 − C_S` ⇒ the shortfall reads as ice ⇒ damps harder.
Small `τ` locks the loop in.

##### Attempt: solid-halo guard — 10× better, still 3.6× short (job 20429954)

Brinkman disabled where `C_S > 1e-3` (≈ `R + 4.2` cells, sized from the `C_S` tanh to
cover the measured `R+3` contamination; the geometric cutoff `1e-12` would reach `R + 26`
cells and kill ice rigidity domain-wide), guarded by `φ_s < 0.9` so genuine ice beside the
grain stays rigid. `U_t`: −0.028 → **−0.279**. Stage-A gate bit-identical.

##### ⛔ ROOT CAUSE FOUND — the grain leaves a PHANTOM-ICE WAKE

Measured in an **ice-free** deck, outside the halo, with the grain at `y = 3.42`:

| region | mean `φ_s` | max |
|---|---:|---:|
| 0–1 d **above** (wake) | 0.0145 | 0.967 |
| 1–3 d above | 0.0308 | 0.993 |
| 3–6 d above | 0.0212 | 0.990 |
| 0–1 d **below** (undisturbed) | **0.0007** | 0.057 |

**Cells the grain vacates are not refilled to `C_L = 1 − C_S` fast enough and read as
ice.** The asymmetry — everything above, nothing below — is conclusive. **97 % of the
domain carried `φ_s > 0` at a mean of 0.022**, i.e. ~19 % spurious drag *everywhere* at
`τ = 1e-3`, on top of concentrated wake cells. No `τ` and no halo radius can fix a
distributed drag. This is the same class as the known moving-IBM wake mass-conservation
artifact in this codebase.

##### Attempt: Yang-style hard threshold — ✅ ice solved, ❌ particle now DECELERATES (jobs 20445849–20445851)

`VOF_DIFFUSE_ICE_PENAL_THRESHOLD` (flag-gated, **undefined by default** so Stage A / Yang /
melting-RB keep the graded mask they were validated with) applies the reference method's
own rule — §G.8: Yang/AFiD use *"η = dt penalty **plus hard `u = 0` for `φ_s > 0.9`**"*.
Damped set drops **97 % → 0.156 %** of cells.

| | result |
|---|---|
| Stage-A gate | ✅ `λ = 0.2184495545645199`, enthalpy `4.975384548799866e-10` — **bit-identical, 8th consecutive** |
| **ice rigidity @ `τ = 1e-3`** | ✅ **0.060 %** leakage (vs 14.06 % at `τ = 0.391` graded) — **SOLVED** |
| particle | ❌ **decelerates monotonically** |

`U_y` vs time: `t=2: −0.833` → `t=5: −0.215` → `t=8: −0.110` → `t=12: −0.071` →
`t=15: −0.052`. It **starts near the correct ~1.0 and then strangles**; `sd = 0.097` on the
mid-fall mean confirms it never reaches a terminal velocity at all.

**Hypothesis for the next session (NOT yet tested):** the threshold makes wake cells that
cross 0.9 *fully* rigid (`φ = 1`) rather than partially damped, so the phantom-ice wake
becomes a **rigid plug above the descending grain**. In an incompressible box the return
flow must pass through that region; as the wake lengthens the obstruction grows, which
matches the monotonic deceleration. The graded mask + halo guard instead gave a *steady*
−0.279 — a different failure, consistent with distributed rather than blocking drag.

##### Historical state of B.2 at handover — superseded by the 2026-09-07 closure audit

| check | status |
|---|---|
| 1 stationary while locked | ✅ exact, 2-D and 3-D |
| 2 impulse-free release | ✅ 0.24–0.68× `a_free` |
| 3 meltwater tracer | ✅ 3.25e-3 in 3-D, bounded |
| 5 budgets | ✅ `net = 0`, enthalpy 2.7e-4, salt 0 |
| **4 settle + wall collision** | ⛔ **OPEN** — works in the *standalone* gate (`U_t = −1.0005`, clean landing at `y = 0.4996`) but not with a valid `darcy_tau` |
| ice rigidity | ✅ `τ = 1e-3` + threshold ⇒ 0.060 % |
| Stage A untouched | ✅ 8 consecutive bit-identical gates |

**The single open item is the phantom-ice wake.** Everything else in Stage B is validated.
The fix must make the ice indicator independent of the CH refill lag behind a moving
solid — e.g. refill `C_L → 1 − C_S` in cells the solid vacates (a solid-motion source for
`C_L`, the direct analogue of the existing `VOF_DIFFUSE_BAND_MASS_RESTORE` contact-line
work), rather than any further tuning of `τ` or of the mask geometry. **Three successive
attempts to fix this by masking or by `τ` have each revealed the next layer; the wake
refill is the first candidate that addresses the cause rather than the symptom.**

##### Historical deck recommendations — withdrawn

**Do not apply the tau values in this historical table.** The validated 3-D
configuration uses actual-ice transport plus the hard threshold at `tau=1e-3`.
Legacy 2-D arms are outside the new flag's supported configuration. A real
scenario deck requires scenario selection and approval; no production deck was
run or EOS fitted in the closure audit.

| deck | change |
|---|---|
| `StageB_release3d` | `darcy_tau` `1e-4 → 0.391` |
| `B14_*`, `Settle_validation` | `darcy_tau` `1e-4 → 0.0076` (2-D) / `0.391` (3-D) |
| **§B.3 production block below** | `darcy_tau` `1e-4 → 0.174`, plus the four errors already listed (`Cn`, `Pe` tracer, `conc_init_type`, `F_release`) and the no-flux thermal BC |

> **Scope of the consequence.** Every settling number this project has ever produced is
> affected — B.1.3-fn's "settles after release", the B.1.4 arms, and B.2. None of them are
> wrong *about what they tested* (lock, release timing, impulse-free handoff, budgets,
> conservation — all of which remain valid), but **no settling velocity or wake from this
> code is currently trustworthy**, and §B.3's science output is exactly that.

> **Cost note for §B.3.** Two full 250 t.u. shakeout runs at 64 ranks cost ~28 h wall
> (~1800 core-h) before check (2) was even reachable, entirely because of a boundary
> condition. The B.3 production deck carries the **same** no-flux thermal BC — at
> production grid (35.4 M cells) that mistake would cost ~100× more. Fix it there before
> submitting anything.

```ini
# DEBUG: short fat column, one particle, coarse but 3D, reduced Sc.
# d = particle diameter; here d/Dx = 16 (debug-only; production uses >=24).
[geometry]
xmin = 0.0
xmax = 4.0              # Lx = 4 d
ymin = 0.0
ymax = 10.0             # Ly = 10 d (vertical)
zmin = 0.0
zmax = 4.0              # Lz = 4 d

[grid]
NXM = 64
NYM = 160
NZM = 64                # ~0.65 M cells

[simulation]
time_max = 20.0
default_dt = 1.0e-3
cfl = 0.3

[flow]
Re = 100.0              # Ga-scale settling

[vof]
init_type = 28

[conc]
NConc = 3
Pe = {700, 700, 0}      # Sc=7 debug; tracer non-diffusive
conc_init_type = {32, 33, 34}

[particle]
F_release = 0.5
# all Stage-A + Stage-B flags on.
```

**What to check:** (1) particle stays put (velocity < 1e-6) while in ice with `darcy_tau` small and
melting off; (2) clean, impulse-free release when the melt front passes (no velocity spike); (3)
meltwater tracer is injected only at the interface and is globally conserved; (4) collision with the
bottom wall fires without blow-up; (5) heat/salt/tracer budgets close. Fix bugs here; do not treat
this coarse low-Sc run as a result.

### B.3 Full ECCO production inputs — scenario design and two precursors

**Design rewritten 2026-09-11; results updated 2026-09-27. Bounded P1 feasibility
job 20906327 COMPLETED; full P1 remains OPEN.** The records below
supersedes unresolved selections and linear-EOS starting defaults in the design.
Stage A and B.2 are closed with the limitations recorded in the
[Stage-A summary](YANG_VALIDATION_SUMMARY.md) and
[B.2 closure report](YANG_ECCO_B2_CLOSURE_REPORT.md). This section replaces the old
B.3 in full. It controls the forward plan wherever the historical critical-path,
resource, or input recommendations elsewhere in this roadmap disagree with it.
The proposed sequence is **two cheaper simulations, then one production case**;
there is no automatic parameter sweep or 48-cells-per-diameter run.

#### B.3 P1b redesign — 2026-09-28 (controlling)

**Full record: [YANG_ECCO_B3_P1b_SHELF.md](YANG_ECCO_B3_P1b_SHELF.md).**
Inputs are in `ECCO_TESTS/StageB3_scenarios/P1b_shelf/`; output goes to
`/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/`. This subsection
controls wherever B.3.1–B.3.6 below (the 2026-09-11 design) disagree with it.

**Why P1a is superseded.** Three findings drive the redesign.

1. P1a's net freezing is a model artifact, not the physics of its inputs:
   - The exact binary-Stefan solution for those inputs melts (λ = 0.0112, interface
     θ = 0.046). A 1-D enthalpy model of the finite column also melts.
   - With `liquid_referenced_salinity = 0` the liquidus reads F·s. That error is about
     0.5 × `liquidus_slope`: invisible at Yang's 0.014, but θ ≈ 0.5 at P1a's 1.01.
   - So the band pinned at θ ≈ 0.44, above the θ = 0.2 threshold where the thin cold
     ice outdraws the water. `F∇s` also piled salt to s/F = 3.2 on the ice side.
2. Waiting for a trapped grain to melt out is unaffordable:
   - The front must rise 0.709 d (shell trigger) before release.
   - At d = 1 mm that takes 3–6×10⁴ code units (0.5–1 h) for Antarctic shelf water,
     and 1.2–1.5×10⁴ for Greenland Atlantic Water.
   - That is about 0.2–1×10⁶ SU at P1a's timestep, and still 2–8×10⁴ SU at the
     explicit Cahn–Hilliard limit dt ≈ 0.021. Yang production ran at 96% of that
     same limit.
   - Only the grain size helps: cost scales as d^1.5.
3. Release timing is ill-conditioned (B.2) and not a deliverable anyway.

P1a's rigidity "failure" was a small-denominator ratio in a still box: the domain
maximum was 2 µm/s.

**P1b scenario (user decisions 2026-09-28).** The release is the initial condition:
- **Grain:** a free (`F_release = -1`) quartz grain, d = 0.2 mm, starting 0.35 d
  below the ice interface of an Antarctic ice-shelf base.
- **Water:** 700 dbar, mixed layer CT = 0 °C, SA = 34.85 g/kg, u* = 0.5 cm/s.
- **Interface state** (three-equation model, Γ_T = 0.011, Γ_S = 3.1×10⁻⁴):
  21.30 g/kg, −1.670 °C, 35 m/yr.
- **Sublayers:** δ_T = 12.7 d, δ_S = 3.2 d, initialized as erf profiles (new Conc
  init 36/37/38) with the tracer equal to the meltwater fraction.
- **Ice:** 2 d, isothermal at the interface temperature, insulated top.
- **Melting ON:** with the band in equilibrium at t = 0. The band treatment was decided
  by two pre-flight rounds against the exact solution at P1b's parameters
  (λ = 0.019382):
  - **Round 1 (job 20950025):** the default operator ∇·(F∇s) piles salt into the band.
    With liquid referencing λ is +94%; without it, −9% with the wrong interface state.
    Both fail.
  - **Round 2 (job 20971700):** the **Yang operator** (`yang_salt_transport = 1`, salt
    field = liquid salinity) gives λ within 0.2% and the correct interface state. Its
    ∫(F+δ)S drift is 1.8×10⁻⁵, which misses the pre-declared 10⁻⁸. It was adopted as a
    recorded deviation; see the dedicated file §8.2.
  - **Separate finding:** `liquid_referenced_salinity = 1` divides by C_L inside resolved
    grains, where C_L → 0 (s/F ≈ 1,564). The buoyancy there is then about 10³× too large,
    and the pressure loop stalls. Never combine it with `LAG_PARTICLE_RESOLVED`.
- **Scaling:** Re = 0.812, Pe_T = 10.43, Pe_S = 56.81, St = 0.0200, G* = 36.78.
- **EOS:** `EOS_NONLINEAR` refitted at 700 dbar. The interface water sits near its
  density maximum, so a linear EOS is inadequate here.
- **Box:** 8 × 22 × 8 d at 24/d (19.46 M cells), T_END = 30 (0.82 s).
- **Control:** a matched ice-only control isolates what the grain does to the layer.

**Cost:** about 1.9k SU each for P1b and the control (2.5k each with 30% margin).
P2 at 32/d is about 6k SU after analysis.

**Result 2026-10-01:** the segment is complete. Highlights are in the header block; full
results are in `StageB3_scenarios/P1b_shelf/RESULTS_20977443.md`.

**Submitted 2026-09-30:**
- P1b-Yang **20973340** failed at t = 0.738 (Yang salt field unstable at the moving grain).
  It was resubmitted with the default operator as **20977443**, with control **20977444**.
  The Yang control 20973341 is kept as an operator reference.
- Original submission: P1b **20973340** and control **20973341**;
- each 1 node, 128 ranks, `shared`, 24 h;
- binary `build_v2`;
- pre-flights cost about 62 SU in total.

#### B.3 P1a result — 2026-09-27 audit (superseded by P1b on 2026-09-28)

**The bounded feasibility segment completed normally; full P1 is not accepted.**
Job **20906327** finished with exit `0:0` on 2026-09-25, 15:44:10–17:09:50 local
Anvil time. Cost is **182.75556 SU = 128 CPUs × 5140 s / 3600**; do not add Slurm
step charges to the main job row. This is **47.59% of the 384-SU cap**. Current
`mybalance`: **71,516.5 SU remaining**. The cancelled original job used zero CPUs.

The reproducible [result report](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/RESULTS_20906327.md)
and [compact analysis](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/results_20906327/summary.json)
include all-field audits, time-series CSVs, accounting, figures and provenance.
Raw output remains under `ECCO_StageB3/P1_Attempts/B3_P1_20260925_20906327/`; the complete
analysis is under `ECCO_StageB3/P1_Attempts/B3_P1_analysis_20906327/`, both beneath
`/anvil/scratch/x-mjalabert/`. No raw fields were removed or copied into home.

| Quantity / gate | Result |
|---|---|
| Completed duration | 5,032 steps; t=10.0036050786 = 0.619943815 s |
| Stored fields | 401 Data, 401 Resume, 401 Particle frames; 284.64 GiB apparent raw size (288G allocated by `du`), above the 250-GiB plan |
| Motion / release | Exactly zero displacement and speed; y=8.5. First computed shell fraction 0.0536644 → 0.116676, below 0.709; no release |
| Phase evolution | Ice 31.4745276870 → 31.4808039071 d³, **+0.00627622009 d³ (+0.0199406%)**; monotonically increasing at all 805 profile times |
| Salt budget | **Pass:** max relative 2.04281e−14; absolute 2.61480e−12 |
| Signed tracer/net-phase budget | **Pass:** max relative 4.34326e−5; absolute 1.11441e−7 |
| Boundary-corrected heat | **Pass:** max relative 1.26602e−10; absolute 1.63864e−7; final error / accumulated boundary heat = 1.85308e−5 |
| Numerical field health | All 401 saved phase/velocity/scalar frames finite; EOS range covered; logged unplaced residual ≤4.46e−19 |
| Ice rigidity | **Relative criterion fails:** all 400 nonzero-time fields exceed 1%; worst **3.113915%** at t=2.599605, final **1.221667%**. Maximum absolute bulk-ice speed **7.26e−8 m/s** |
| Post-release / refined comparison | Handoff, escape, settling, wake and P2 sensitivity **unevaluated** |

**New physical/model finding:** local melting and freezing coexist, but the
net tendency over this short startup is **freezing**, not the melt-out assumed
by a simple heat-only forecast. Saved-time source integrals are approximately
+0.061998 and −0.068271 d³; these are sampled diagnostics, not exact RK totals.
The first positive-time shell value is the first computed one; the t=0 zero
in `release.dat` is an initialization placeholder. Do not extrapolate the shell
fraction or this short net-freezing trend into a release clock or a claim that
later melt-out is impossible.

The source adds signed `m` to tracer, including negative values for freezing.
The physical tracer range is **−0.0176284 to +0.0113315** outside sediment
support. No saved bulk-liquid cell (`C_L>0.95`, `C_S<0.05`) is below −1e−10;
negative values occur in the interface/ice region. Thus its signed conservation
check passes, but it cannot everywhere be interpreted as a nonnegative
meltwater-origin fraction. Resolve that interpretation before entrainment/mixing
claims; simply clipping negative tracer would destroy the measured balance.

The rigidity comparison uses the recorded B.2 mask (actual ice >0.95,
`C_S<0.001`). Absolute ice speed is tiny, but the persistent 1.18–3.11% relative
result does not meet the declared 1% bar. Diagnose it before continuing; no
Darcy retuning or silent acceptance-threshold change was made. Ice integrals
exclude `C_S>=0.05`, matching the solver; independent same-time field/profile
ice totals agree to **5.33e−14 d³**. C_S is a shifted support indicator, so its
raw integral must not be interpreted as the geometric sphere volume.

**Measured planning update:** throughput is **16,916.9 cell-steps/core-second**,
or **18.269 SU per code time unit** under this output policy. Applying it to
the old one-diameter conduction scenario costs about **665,627 SU** before
contingency; this is still a conditional planning scenario, not a measured
melt-out cost. A matched *short* P2 would be about **866 SU / 1,126 SU with 30%
contingency** if throughput transfers; it would not provide the required
post-escape comparison. Full-field output during long locked intervals needs
replanning: it took ~16% of solver wall time and produced 284.64 GiB here.
The measured wall time plus 30% margin is **111.37 minutes**, so an identical
bounded repeat could plausibly fit `debug`'s 2-hour limit after revising its
stop timer; this says nothing about long P1/P2 feasibility. No repeat is submitted.

**Continuation decision:** preserve the run; investigate signed freezing/tracer
semantics and the rigidity miss before spending on a long continuation or P2.
The matching restart is `Resume_400.h5` (transport version **2**, step 5030,
t=**9.99960507859**), two timesteps before the final profile; all matching data
and histories are retained. This audit verified files/version, not a new restart
run. No simulation was automatically launched because the grain failed to release.

#### B.3 submission record — 2026-09-25 (historical plan and provenance)

**Authorization and scope.** The user requested proceeding with B.3/P1 and
accepted the proposed idealized physical case. This satisfies the earlier P1
approval boundary. It does not authorize P2 or full production. The allocation
query immediately before submission returned **71,699.3 SU**, with no active
jobs. The old 77,010-SU balance and synthetic timestep advice below are historical.

Reproducible inputs, EOS fit, physical/scaling record, job script and audit are
in [`B3_P1`](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/README.md), particularly
[`case.json`](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/case.json) and
[`prepare.py`](../../../PARTIES/testcases/ECCO_TESTS/StageB3_scenarios/P1a_feasibility/prepare.py).
No shared solver source or existing executable was changed.

| Selected physical quantity | Value |
|---|---:|
| Grain diameter / mineral density | 0.001 m / 2650 kg/m³ |
| Initial water / maintained bottom Conservative Temperature | 0 °C |
| Interface Absolute Salinity / stable gradient | 34 g/kg / 0.01 (g/kg)/m |
| Fixed reference sea pressure | 0 dbar |
| Initial ice / top temperature | CT freezing at SA=34: −1.84876067498 °C |
| Reference density | 1027.17207243 kg/m³ |
| Fixed model ν / κT | 1.8e−6 / 1.4e−7 m²/s |
| CT heat capacity / constant latent heat | 3991.86795712 / 333500 J/(kg K), J/kg respectively |
| Sc / Pr / Le | 70 / 12.85714286 / 5.44444444 |
| Modeled salt diffusivity / physical comparison constant | 2.57142857e−8 / 1e−9 m²/s |
| Uref / tref | 0.0161363092 m/s / 0.0619720401 s |
| Re / Ga / G* | 8.96461621 / 69.16344762 / 37.67563614 |
| PeT / PeS=PeC / St | 115.2593513 / 627.5231350 / 0.0221289610 |

These material constants are declared idealizations, not measured site values.
The equal ice/water thermal diffusivity and unheated rock approximations remain.
The ambient salinity change over the 8-mm initial liquid depth is only
0.00008 g/kg; the chosen gradient is preserved, not artificially amplified to
make millimetre-scale ambient stratification strong. The fitted local ambient
`N²=7.5400e−5 s−2` is stable. Meltwater supplies much stronger local density changes.

**EOS gate changed the design.** The full-dilution linear density fit has
`dρ/dCT=+0.0186151 kg/(m³ K)`: it reverses the thermal buoyancy direction in
ambient seawater. It is rejected despite small total-density error. A
free density-only nonlinear fit also masks large errors in thermal response.
Use the **existing `EOS_NONLINEAR`, q=2** instead: pin the density-maximum
locus to GSW, fit curvature to the thermal derivative, then fit the haline
coefficient and additive constant. The constant is absorbed by pressure.
GSW 3.6.23 exact Gibbs density, after explicit CT-to-in-situ-temperature
conversion, covers SA=0–34.00022 g/kg and CT=−1.84876–0.03 °C; a separate denser
grid checks the fit. This includes diffuse-band/supercooled states as a numerical
coverage check, not a claim that every sampled liquid state is thermodynamically stable.
The maximum errors are **0.0698261 kg/m³**, **0.256115% of the haline density
scale**, **0.00389982 kg/(m³ K)** in thermal derivative, and
**0.0131292 kg/m³ per (g/kg)** in salinity derivative. Over SA≥28 the maximum
relative thermal-derivative error is **6.402%**. The density-maximum locus error
is at most 0.03913 K; signs immediately around that locus are approximate.

The selected input coefficients are `eos_betaT=0.0009466648317134`,
`eos_betaS=1`, `eos_Tmd0=3.268946802839`, `eos_Tmd_slope=−4.207161725279`.
Scaling beta is the **explicit** haline coefficient divided by rho0; the
quadratic contributes an additional local salinity derivative. `richardson`
is inert in this build. Thus the linear-EOS rows and example in B.3.2/B.3.4
remain the original design starting point, **not the selected P1 configuration**.
The checked approximation is suitable for this bounded feasibility observation;
its errors remain part of the physical interpretation, not a TEOS-exact claim.
Any matched P2/production must retain the same selection.

The freezing secant is exact at SA=0 and 34, with freshwater CT freezing
0.0179473461 °C: `T_melt=liquidus_slope=1.009707771432` and therefore
`theta_liquidus(1)=0`. Maximum interior error is **0.020810 K** across dilution.
Sources: [GSW density](https://teos-10.github.io/GSW-Python/density.html) and
[GSW freezing/temperature functions](https://teos-10.github.io/GSW-Python/gsw_flat.html).

**Affordability finding.** The selected case has much smaller St and larger PeT
than the synthetic B.2 case. The initial liquid sensible heat melts only
**0.1770d** of a planar layer. A simple heat estimate for one-diameter retreat,
using that reservoir first and an 8.5d mean conduction path, requires about
**36,435 nondimensional units / 2,258 physical seconds**. At dt=0.002 the two
historical throughput anchors imply **357,597–2,417,449 SU**, before contingency.
This is a planning estimate, **not a release-time prediction or rigorous bound**:
cold-top extraction reduces heat supply; convection can increase it; dilution
changes freezing. It does not justify spending the account on a full melt-out
attempt, or increasing St/temperature just to obtain release.

**Submitted observation.** Job **20906327**, `ECCO_P1_feas`, is a bounded first
P1 segment on the full planned `96×240×96`, `4×10×4` grid with the original
sphere at `(2,8.5,2)`, radius 0.5, interface y=8. It uses `max_dt=0.002`,
`default_dt=0.0002`, `time_max=10` (0.61972 s), fields every 0.025, profiles
every ten steps, and 50 increment-ramp steps (target duration 0.1 only if the
cap is reached). The source relaxation scale is 0.90046; speed bound 5 gives
advective CFL 0.24. Adaptive CFL and particle/contact limits remain active.

At the cap, 5,000 steps cost 98–663 SU, or **128–863 SU with 30% contingency**.
The job has an intentionally smaller **hard cap of 384 SU** (128 cores × 3 h),
with normal checkpoint stop requested at 170 minutes. Reaching t=10 is therefore
conditional on throughput. Storage allowance is about 250 GiB. This segment
measures actual startup, heat/melt direction, conservation, geometry and cost;
it does **not** promise a post-escape interval or complete P1 acceptance. No
automatic long continuation is authorized by failure to release. A matched
t=10 P2 costs approximately 465–3145 SU (605–4089 with contingency), but would
likewise not establish the required post-escape comparison.

The isolated ordinary build is
`/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_build_20260925`; source hashes and
binary hash are in `manifest.json`. Its header is identical to
`B2_closure/Boundary.validation.h`, retaining the nonlinear EOS and fixed release
ramp. Shell syntax, Python compilation, physical sign checks, independent EOS
holdout evaluation, liquidus endpoint consistency, and build completion passed.
The new output audit also ran successfully on archived B.2 job 20492853; that
checks the audit reader, not P1 physics. Its smoke output is archived beside the build.
Runtime inputs, executable, manifest, audit and case record are archived to
`/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_20260925_20906327` when the job starts.
The initial scheduler state was **PENDING**; no runtime acceptance is claimed.

**Storage/queue update (user instruction, 2026-09-25).** All future B.3
simulation output, builds, staged submissions and scheduler logs use
`/anvil/scratch/x-mjalabert/ECCO_StageB3/`. Historical B.2 archives remain in
`ECCO_StageB`. Original job **20906230** was verified PENDING, held, and cancelled
before execution; replacement **20906327** was submitted from
`ECCO_StageB3/P1_Attempts/B3_P1_submit_20260925`. The isolated build was moved and its binary
SHA-256 reverified unchanged. The live Slurm configuration lists no `short`
partition; the short CPU option is `debug`, whose `part-debug` QoS allows at most
**2 hours**, 2 nodes / 256 CPUs, and one job per user. It cannot accommodate this
3-hour job as configured, so the replacement retains **shared**, 128 cores,
3 hours, and the **384 SU cap**, with identical physics and inputs.

#### B.3.1 Scientific scope and physical case

> **2026-09-28:** B.3.1–B.3.6 are the 2026-09-11 design, kept as history and as the source of
> the unchanged rules (scaling form, budgets, output policy, cost method). Where they conflict
> with **§B.3 P1b redesign** above, P1b controls. The conflicting items are the grain size (now
> 0.2 mm), release as a melt-out outcome (now the initial condition), the p = 0 and
> linear-EOS defaults (now 700 dbar with `EOS_NONLINEAR`), the cold-top/warm-bottom Dirichlet
> walls (now insulated), and the geometry and timestep tables.

The selected setting (user decision, 2026-09-11) is an **idealized local column of
stratified water beneath stationary, initially fresh ice, containing one resolved
sediment sphere**.
It retains the sediment-from-ice objective. It is not yet a selected sea-ice site
or ice-shelf site, and does not represent brine drainage, ice mechanics, a moving
ice body, or a regional ocean circulation. The first campaign measures the
coupled meltwater field, particle settling, wake evolution, and scalar transport.

**Release is a modeled initialization process, not a quantitative clock.** B.2
did not establish converged escape timing. Do not report absolute release time,
escape duration, or release-to-contact time as physical predictions; increasing
the particle count does not remove the bias. Compare post-escape particle speeds
at matched heights, and inspect the surrounding temperature/salinity/tracer
state as well: shifting the time origin alone cannot correct a different melt
history. If the scientific question requires a prescribed, reproducible release
state, the alternative is an explicitly **already-free grain** with
`F_release=-1`, fully outside the ice and its diffuse support. That is a separate
initial condition to choose before the runs, not a silent substitution for melt-out.

The physical case must supply the following quantities before numerical input
files are finalized. These are explicit unresolved selections, not values copied
from the synthetic B.2 deck:

| Physical selection | Required definition / use |
|---|---|
| Grain | Diameter `d` in metres and material density `rho_p`; enter `rho_s=rho_p/rho0`, rather than treating the old value 2.65 as a universal density ratio |
| Ocean temperature | Initial water temperature and maintained bottom temperature `T_b`; the baseline uses the same value for both |
| Salinity | Absolute Salinity immediately below the initial ice interface `S_i`, and stable vertical gradient `gamma_S=-dS_A/dy_dim >= 0` |
| Depth | Local reference sea pressure `p0`; use a fixed-pressure EOS for this small column, with the pressure chosen for the site |
| Ice temperature | Initial/top-wall temperature `T_i`; the initial baseline uses ice at the ambient interface freezing temperature |
| Material properties | Consistent `rho0`, liquid `nu`, thermal diffusivity `kappa_T`, heat capacity `cp`, latent heat `L`, and the physical salt diffusivity for comparison with the modeled one |
| Scope of the approximation | Initially fresh ice, Boussinesq equal ice/water reference density, constant viscosity, and the declared ice/water thermal-diffusivity ratio |

Use one selected physical case for both precursors and production. Do not adjust
`St`, salt forcing, gravity, or water temperature simply to make the small case
release sooner. If its physical melt time is unaffordable, revise the scientific
initial condition or scope explicitly before committing more time.

#### B.3.2 One consistent scaling, EOS, and liquidus

**Length unit is the grain diameter `d`; velocity uses a haline buoyancy scale.**
This follows the polar-regime finding of B.1.1 without mixing the Yang thermal
free-fall units with the B.2 synthetic `g=d=1` convention. Define:

```text
S_ref = S_i > 0
T_ref = T_freeze(S_i, p0)
DeltaT = T_b - T_ref > 0
theta = (T - T_ref)/DeltaT;       s = S_A/S_ref

delta_h = beta * S_ref            # beta in (g/kg)^-1, from the selected EOS
U_ref = sqrt(g_dim * delta_h * d)
t_ref = d/U_ref
G_star = g_dim*d/U_ref^2 = 1/delta_h

Re = U_ref*d/nu
Pr = nu/kappa_T
Pe_T = Re*Pr
Pe_S = Re*Sc_sim
Pe_C = Pe_S                       # first model: tracer and salt have equal diffusivity
St = cp*DeltaT/L
rho_s = rho_p/rho0
Ga = Re*sqrt((rho_s - 1)*G_star)   # Ga is NOT equal to the input Re
```

The particle input is `grav={0,-G_star,0}`. In the scalar buoyancy assembly,
`Velocity.c` normalizes `grav` to a direction and adds the scalar density anomaly
along it. Therefore the **linear** path must use

```text
richardson = {-G_star*alpha*DeltaT, G_star*beta*S_ref, 0}
            = {-alpha*DeltaT/(beta*S_ref), 1, 0}
```

where `alpha` and `beta` belong to the same density fit and reference state.
Increasing temperature then accelerates water upward, increasing salinity
accelerates it downward, and the tracer is passive. The positive-temperature /
negative-salinity coefficients printed by `eos_polar_refit.py` use the opposite
buoyancy convention and **must not be pasted directly into `richardson`**.
This mapping was checked against the momentum assembly, not inferred from the
old comments on the input keys.

**EOS choice:** start from the existing B.1.1a recommendation, with
`EOS_NONLINEAR` **off**, and fit a linear density model including an additive
constant. Check its density-anomaly error and buoyancy directions over the
actual liquid range, including dilution toward `S_A=0`, at the selected `p0`.
The fit constant sets the reference pressure/density convention; it is not a
fourth active scalar. Report error relative to `rho0*delta_h`, along with local
errors where small buoyancy differences control the science. The July ambient
polar fit is useful prior work, not a completed fit for this case. If a linear
model is inadequate, select and check the alternative before P1; do not switch
EOS between the precursor pair and production.

Use [TEOS-10/GSW density routines](https://teos-10.github.io/GSW-Python/density.html)
and the corresponding
[freezing-temperature function](https://teos-10.github.io/GSW-Python/gsw_flat.html#gsw.CT_freezing)
with Absolute Salinity, Conservative Temperature, and sea pressure consistently.
Document the temperature convention when assigning heat capacity and latent heat.
For the solver's linear liquidus, construct

```text
T_freeze(S_A, p0) ~= a_f - m_f*S_A,      m_f > 0
T_melt = (a_f - T_ref)/DeltaT
liquidus_slope = m_f*S_ref/DeltaT
theta_liquidus(s) = T_melt - liquidus_slope*s
```

Constrain/check the fit at `S_i` so that `theta_liquidus(1)=0` within the stated
fit error, and report its error across the dilution interval. At nonzero `p0`,
use the pressure-dependent freshwater freezing point when choosing `a_f`.
**`T_melt=0` is generally wrong with this temperature origin.** Neither the Yang
value `liquidus_slope=0.014` nor the old B.3 value 0.020 is a polar default.

**Diffusivity model:** the first proposed campaign uses **`Sc_sim=70`**, held
fixed across P1, P2, and production, with `Pr` obtained from the chosen material
properties. This is a reduced salt-Schmidt-number model, subject to resolution
and cost checks, not the physical ocean diffusivity. Check that the chosen
`Sc_sim > Pr` if the intended regime requires salt to diffuse more slowly than
heat. `Le=Pe_S/Pe_T=Sc_sim/Pr`; it is not equal to `Sc_sim`. Do not infer an
oceanic flux law by extrapolating two Schmidt numbers without evidence. If this
model is unaffordable or unresolved, revise it explicitly; reducing `Sc` only
in the precursor would stop that run from testing the proposed production physics.

#### B.3.3 Geometry and the two cheaper simulations

All cases use the supported **3-D uniform Cartesian grid**, gravity in `-y`,
periodic `x,z`, and no-slip top/bottom walls. The initial ice thickness is `2d`;
the sphere has radius `0.5d`, is centered horizontally, and its lower surface
initially touches the nominal planar ice interface. It is initially locked.
Check its initial shell fraction and actual ice geometry from the first output.

| Input / purpose | P1: scenario precursor | P2: refined precursor | First production candidate |
|---|---:|---:|---:|
| Domain `Lx x Ly x Lz`, in `d` | `4 x 10 x 4` | Same as P1 | `8 x 24 x 8` |
| Cells per diameter | 24 | 32 | 24, conditional on P2 |
| `NXM, NYM, NZM` | `96, 240, 96` | `128, 320, 128` | `192, 576, 192` |
| Physical cells | 2,211,840 | 5,242,880 | 21,233,664 |
| Initial interface `y_int/d` | 8 | 8 | 22 |
| `vof_slab_x0` | 0.8 | 0.8 | `22/24 = 0.9166666666666667` |
| Sphere `(x,y,z,R)` | `(2,8.5,2,0.5)` | Same as P1 | `(4,22.5,4,0.5)` |
| `Cn=0.75*dx` | 0.03125 | 0.0234375 | 0.03125 |
| `melt_band_eps=2*Cn` | 0.0625 | 0.046875 | 0.0625 |
| `Pe_CH=0.9/Cn` | 28.8 | 38.4 | 28.8 |
| Maximum timestep | `DT_CAP`, selected below | `DT_CAP/2` | At most the accepted P2 cap, rechecked on the larger domain |

**P1** tests the actual chosen salinity forcing, liquidus, heat supply, coupled
transport, and cost in a small box. It should reach a useful post-escape interval
if the physical melt time fits its budget. **P2** repeats that same dimensional
case and geometry with finer space and time resolution. Use a fresh initialization
on its own grid; do not interpolate a P1 checkpoint and call it a matched run.
The pair is a **combined resolution-sensitivity check**, not separate estimates
of spatial and temporal convergence order. The `Cn`, mobility, and melt-band
changes are part of the declared interface-refinement path.

The choice `melt_band_eps/Cn=2` retains the ratio of the accepted **coupled** B.2
deck. The `Pe_CH=0.9/Cn` refinement rule is an explicit extension of its coarse
value, not something the ice-free settling controls validated. Check melt volume
and interface/scalar profiles under this refinement. The initial tanh profile has
5–95% thickness approximately `8.33*Cn`, hence **6.25 cells**, not the old claim
of a 3–4-cell band. Resolving the grain does not establish salt or interface
accuracy: P2 must inspect those fields independently.

P1 and P2 have the same `S_i` and dimensional `gamma_S`. Preserve that gradient
in production as well, rather than preserving the same top/bottom salinities
over a different height. The implemented initializer is

```text
s_water(y) = cbd5 + (cbd2-cbd5)*y/Ly      # y and Ly are nondimensional
s_initial(y) = s_water(y)*F_layer(y)

k_s = gamma_S*d/S_ref
cbd5 = 1 + k_s*y_int                    # water value at the bottom
cbd2 = cbd5 - k_s*Ly                    # extrapolated value at the TOP WALL
```

Thus `cbd2` is **not** the salinity immediately under the ice. Require nonnegative
initial salinity and include the deeper production bottom salinity in the EOS
fit range. With initially uniform water temperature the ambient stability is
`N^2 = g_dim*beta*gamma_S`; in code units `N_star^2=G_star*beta*gamma_S*d`.

The production size is a **candidate**, chosen to increase lateral separation
and usable fall distance without immediately returning to the old 35.4-million-
cell / 283-million-cell plan. It is still a periodic array with a bottom wall.
The small-box pair does not certify absence of confinement in the larger box.
Report that geometry and use an observation window clear of the ice and bottom;
claim isolated-grain or terminal behavior only if the saved fields actually
support it. A 32-cells-per-diameter production version would have
`256 x 768 x 256 = 50,331,648` cells and requires a revised cost estimate.

#### B.3.4 Solver configuration and full input specification

Start from the promoted **release-ramp-fixed** source and the recorded
[`B2_closure/Boundary.validation.h`](../../../PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/Boundary.validation.h),
using an isolated build. Preserve its supported momentum/pressure, IBM, scalar,
and boundary configuration, with the explicit EOS selection above. The required
features are:

| Keep enabled | Keep disabled |
|---|---|
| `GRID_UNIFORM`, `VOF_DIFFUSE`, `VOF_IBM`, `LAG_PARTICLE_RESOLVED` | `TWOD_CARTESIAN`, `AXISYM_RZ` |
| `CONC`, `BOUSSINESQ`, `CONC_CENTRAL`, `CONC_FULLY_IMPLICIT`, `CONC_VOF_PHASEWEIGHTED` | `EOS_NONLINEAR` for the selected linear model; `IBM_SCALAR`, `VOF_SCALAR` |
| `PHASE_CHANGE`, `ICE_PENALIZATION`, `VOF_DIFFUSE_ICE_PENAL_THRESHOLD` | `SURFACE_TENSION`, `VOF_WETTING`, `VOF_GRAVITY` |
| `VOF_DIFFUSE_SEDIMENT_CH_MASK`, `VOF_DIFFUSE_SHELL_FRACTION_NONSOLID`, `VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT` | `VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE` |
| `ECCO_PROFILES`, `ACTM`, `LUBRICATION_NORMAL` | `PARTICLE_RELEASE` — this is a different mechanism, not the `F_release` motion lock |

Keep the normal owner-coordinate / periodic minimum-distance path. Never use
the manufactured-remap fixture or the smooth-penalization diagnostic executable
for a scenario. Preserve the validated collision inputs and solver tolerances
from `coupled16_rampfix`; do not tune contact parameters to interpret a wall
impact as an exact nonpenetration or sediment-deposition prediction.

The following is the **complete scenario-defining production override** of
[`coupled16_rampfix/parties.inp`](../../../PARTIES/testcases/ECCO_TESTS/StageB2_closure/B2_closure/coupled16_rampfix/parties.inp).
Unlisted settings inherit that recorded deck, including inactive periodic-face
scalar BC entries and solver tolerances. **Uppercase tokens are intentionally
unresolved numbers, not valid parser expressions.** Resolve them using B.3.2 and
the precursor results; materialize one flat, numeric `parties.inp` before a run.
The parser does not implement inheritance, arithmetic, or token substitution.

```ini
[geometry]
xmin = 0.0
xmax = 8.0
ymin = 0.0
ymax = 24.0
zmin = 0.0
zmax = 8.0

[grid]
NXM = 192
NYM = 576
NZM = 192
ImportGridFromFile = 0

[simulation]
time_max = T_END
output_time_interval = DT_FIELD
output_time_interval_2d = DT_FIELD
default_dt = DT_START
max_dt = DT_PROD_CAP
constant_dt = 0
cfl = 0.3
ghost_nodes = 3
resume = 0
freeze_velocity = 0

[flow]
Re = RE_CASE
vel_init_type = 0
dp_dx = 0.0
ubulk_target = 0.0
startup_time = 0.0
startup_init = 0
startup_velocity = {0.0, 0.0, 0.0}

[vof]
We = 1.0e30
rho1 = 1.0
rho2 = 1.0
mu1 = 1.0
mu2 = 1.0
init_type = 28
vof_slab_x0 = 0.9166666666666667
Cn = 0.03125
Pe_CH = 28.8
darcy_tau = 0.001
weno_order = 5
ch_iter_max = 3000
ch_tol = 1.0e-11

[phase_change]
stefan = ST_CASE
T_melt = TM_CASE
liquidus_slope = LAMBDA_CASE
melt_band_eps = 0.0625
meltwater_tracer = 1
yang_salt_transport = 0
liquid_referenced_salinity = 0

[eos]
# Inert with EOS_NONLINEAR off; no synthetic nonlinear forcing is retained.
eos_q = 2.0
eos_betaT = 0.0
eos_betaS = 0.0
eos_Tmd0 = 0.0
eos_Tmd_slope = 0.0

[conc]
NConc = 3
Pe = {PE_T_CASE, PE_S_CASE, PE_S_CASE}
richardson = {RI_T_CASE, 1.0, 0.0}
V_s0 = {0.0, 0.0, 0.0}
kappa_ice_ratio = {1.0, 0.0, 0.0}
conc_init_type = {34, 35, 0}
cbd0 = 1.0
cbd1 = 0.0
theta_ice = 0.0
cbd2 = S_TOP_EXTRAP
cbd5 = S_BOTTOM
# A*dc/dn + B*c = C. Fields: temperature, salinity, meltwater tracer.
# Maintained cold top and warm bottom; salt/tracer have no wall flux.
BC_AN = {0.0, 1.0, 1.0}
BC_BN = {1.0, 0.0, 0.0}
BC_CN = {0.0, 0.0, 0.0}
BC_AS = {0.0, 1.0, 1.0}
BC_BS = {1.0, 0.0, 0.0}
BC_CS = {1.0, 0.0, 0.0}

[particle]
rho_s = RHO_RATIO
rho_prImp = RHO_RATIO
grav = {0.0, NEG_G_STAR, 0.0}
N_forcing_loops = 1
N_heating_loops = 0
N_impheating_loops = 0
F_release = 0.709
release_shell_cells = 2.0
release_ramp_steps = N_RAMP

[output]
output_vfc = 1
post_processing_switch = 0
```

Production `p_mobile.inp` (first line is the particle count; last column is
**radius**, not diameter):

```text
1
4.0 22.5 4.0 0.5
```

`p_fixed.inp` contains `0`. For P1/P2, substitute the geometry, grid, interface,
particle position, and interface parameters from B.3.3, recalculate `cbd2/cbd5`,
and use their planned duration and output/timestep settings. All physical
coefficients and the selected diffusivities remain the same.

| Token | Numerical definition before materialization |
|---|---|
| `RE_CASE`, `ST_CASE`, `RHO_RATIO` | `Re`, `St`, `rho_s` from B.3.2 |
| `PE_T_CASE`, `PE_S_CASE` | `Re*Pr`, `Re*Sc_sim`; both strictly positive |
| `RI_T_CASE`, `NEG_G_STAR` | `-G_star*alpha*DeltaT`, `-G_star` |
| `TM_CASE`, `LAMBDA_CASE` | Liquidus intercept and slope from the checked freezing fit |
| `S_TOP_EXTRAP`, `S_BOTTOM` | `cbd2`, `cbd5` from the geometry-specific gradient mapping |
| `DT_START`, `DT_PROD_CAP`, `N_RAMP` | Selected and documented using B.3.5; actual adaptive ramp duration must be recorded |
| `T_END`, `DT_FIELD` | Observation window and output cadence selected using B.3.5; not the old arbitrary 100 and 5 |

`darcy_tau=1e-3` with the hard threshold is the **starting numerical candidate**
from B.2. Its dimensional damping time is `tau*t_ref`; the changed scaling and
active polar forcing mean rigidity must be rechecked. The rejected 0.391/0.174
values and a `Re*dx^2` tuning rule are not alternatives. Retain `mu2/mu1=1`, as
in the accepted coupled deck; do not reintroduce a viscosity ratio of 100 as an
unmeasured rigidity repair. Likewise `F_release=0.709` and a two-cell shell are
the inherited modeled trigger, not a universal geometric detachment criterion.

The baseline deliberately retains **unit ice/water thermal diffusivity ratio**,
matching the supported B.2 heat-budget audit. This is a declared model
approximation, not an assertion about real ice conductivity. A physical ratio
requires the corresponding boundary-flux audit to be extended and checked before
using it. `N_heating_loops=N_impheating_loops=0` also preserves the unheated-rock
model; the sphere does not melt or generate tracer.

`ECCO_PROFILES` already writes the necessary distributed profiles every ten
steps and at field output. `post_processing_switch=1` is **not** needed to enable
them. Use `C_S<0.05` to exclude sediment-support extension values from physical
scalar plots; preserve the implemented full-physical-mesh scalar integrals for
budget checks, and liquid weights for means/fluxes.

#### B.3.5 Timesteps, duration, and precursor acceptance

**Select timesteps in the new units.** B.2's 0.01 cap and 0.28% ice-free settling
comparison are evidence for that control only. They do not supply a production
timestep for different gravity, EOS, `Pe`, `St`, grid, or threshold activity.
For P1, set a preliminary `DT_CAP` using estimated fluid/particle speeds, the
advective CFL, and the explicit melt-source relaxation scale
`4*Pe_T*Cn*melt_band_eps`; retain a safety margin below that scale. Check the
actual source rate and particle displacement per step, and leave the existing
particle/contact timestep restriction active. Choose `DT_START <= 0.1*DT_CAP`.
P2 halves both values. Production starts no larger than P2's accepted cap, with
the larger-domain CFL and scalar/interface behavior checked again.

Use the fixed ramp implementation, which ramps velocity **increments**. Record
the intended nondimensional ramp duration and choose integer `N_RAMP` to
approximate it from the timesteps actually taken near release. Doubling the
step count when halving a **cap** does not guarantee a matched duration under
adaptive stepping. Measure and disclose the actual duration; do not use an
unmatched ramp to claim timestep convergence. Matching it still does not cure
the escape-timing limitation.

Before selecting `time_max`, estimate the heat required to uncover the grain
and the net heat supply, including conduction to the cold top. Warm-bottom
Dirichlet forcing avoids the previously exhausted insulated reservoir but does
not guarantee a specified melt rate or eventual release for every physical case.
Choose a **bounded** precursor duration covering the predicted melt-out plus a
useful fall interval, or explicitly record that the proposed melt-out case is
too expensive. During that interval sample post-escape motion only where the
grain is separated from actual ice and the bottom; use `2 < y_p < 6` as an
initial small-box analysis window, subject to the saved ice/wake fields.
No automatic long continuation is triggered merely because release has not fired.

Select full-field cadence to resolve the falling grain (initial target:
`DT_FIELD <= 0.25/|U_p|max`, so displacement between fields is at most `0.25d`).
Particle time series and every-ten-step profiles carry the faster diagnostics.
Choose production duration from the observed settling and wake evolution, and
record its planned physical duration `T_END*t_ref`, wall time, and storage cost.
Do not assume the small column reaches a terminal plateau.

The following are **prospective decision criteria**, not claims that P1/P2 have
passed. Both runs must satisfy the conservation/geometry criteria; P2 supplies
the combined resolution comparison.

| Decision criterion | Required evidence |
|---|---|
| Physical setup | Correct initial geometry, locked grain, stable ambient density profile, active haline forcing, positive finite scalar diffusivities, and a checked EOS/liquidus over the values actually sampled |
| Numerical health | No unplaced-ice/remap residual, unintended ice generation behind the grain, solver failure, or persistent unphysical liquid-region scalar oscillation; no tracer source attributable to rock |
| Lock / handoff / rigidity | Exactly stationary before the trigger; inspect the acceleration at handoff relative to the same buoyant-weight/inertial-mass convention as B.2; retain its `<2*a_free` diagnostic bar. Sample bulk-ice speed and require at most 1% of the domain speed where that normalization is meaningful; also report absolute speed |
| Scalar and heat budgets | Salt relative drift `<=1e-8`; tracer/net-melt mismatch `<=1e-3` once net melt is significant; boundary-flux-corrected relative heat error `<=1e-4`. Report absolute errors too and use the recorded B.2 normalization; investigate errors hidden by small denominators or cancellation |
| Post-escape settling | On a common ice-clear, wall-clear height window, P1/P2 maximum speed difference `<=2%` of the sampled peak speed; compare ambient scalar states too. No timing-based acceptance item |
| Scalar/interface and transport sensitivity | Check melt volume before escape at common physical times, interface position/profiles, and gradients in the physical liquid. Target `<=5%` normalized change in selected nonzero liquid-weighted flux profiles and tracer penetration at comparable states; use absolute scales near zero and report sign changes |
| Resolved scales | Examine thermal/salt/tracer gradient thicknesses and, if turbulent, dissipation-based viscous/Batchelor scales. A resolved particle and a successful budget do not certify DNS of the scalar field |
| Cost and boundaries | Measured steps, solver iterations, memory, throughput, I/O, and usable fall window; identify periodic-image and bottom-wall influence rather than fitting it away |

Agreement of this two-run pair is a screening criterion, not proof of asymptotic
convergence: spatial, temporal, and diffuse-interface errors can interact.
If it fails, do not launch production on the premise that a larger box cures it.
Use the observed failure to revise the resolution/model and costed plan. With
only one precursor, P1 can establish operability and cost; it leaves the
scenario-specific resolution comparison open.

For diagnostics, keep `H = integral(T) - V_ice/St` consistent with the implemented
actual-ice accounting, and subtract the **integrated top and bottom heat fluxes**.
Track salt and tracer on the physical mesh, excluding duplicated high planes.
Use liquid-weighted means, covariance fluxes, and time-integrated dissipation.
Evaluate `kappa_eff=-<v'c'>/(d<c>/dy)` only where the denominator is resolved and
non-negligible. BPE must use the selected EOS and changing liquid volume, with
melting/source/boundary contributions accounted for before any mixing-efficiency
claim. Do not reuse the old endpoint-dissipation formula for `eta`.

One coupled production run measures **total coupled transport**. Attributing
an enhancement specifically to the grain requires a matched ice-only background
control; that is additional science scope, not an effect established by the two
resolution precursors. The accepted Stage-A salty-rate residual remains a model
limitation when interpreting absolute melting rates here.

#### B.3.6 Cost, reproducibility, and the decision to proceed

Estimate cost from the corrected, archived runs and then replace that estimate
with P1/P2 measurements. For example, the fixed coupled job **20492853** used
29.33 SU for 393,216 cells and 8,406 steps; the completed 32-cell-per-diameter
settling control **20477407** used 248.36 SU for 5,242,880 cells and 789 steps.
These imply roughly **31,300 and 4,630 cell-steps per core-second**, respectively.
The spread is material. Neither case includes this scenario's active polar salt
forcing and scalar parameters, so these are planning anchors, not performance
bounds or promises. Values come from `B2_closure/campaign_jobs.json` and the
corresponding scratch profiles; the account balance quoted elsewhere is historical.

```text
SU estimate = N_cells * N_steps / (3600 * measured_cell_steps_per_core_second)
wall hours  = SU / allocated_cores
```

At those two throughput anchors, **10,000 steps** cost approximately 200–1,330 SU
for P1, 470–3,150 SU for P2, and 1,900–12,750 SU for the proposed production grid,
before contingency. This is not a fixed-duration comparison: refinement and
the physical melt time can change the step count substantially. Set a duration,
allocation cap, and storage estimate for each run from its numeric deck; include
at least 30% planning contingency and refresh the available allocation before
submission. Do not promise a cheaper run solely from its cell count.

Archive the flat input files, `p_mobile.inp`, `p_fixed.inp`, selected physical
case/scaling/EOS record, boundary header, source/build manifest, executable hash,
and audit scripts for each run. Outputs belong in unique directories under
`/anvil/scratch/x-mjalabert/ECCO_StageB3/`, never in home. Use an ordinary production
build, not a test fixture. A restart must retain transport checkpoint **version 2**
and all matching data/history files on the same grid and physics; it is not a
mechanism for converting a precursor into the larger-domain production case.

**Proceed in this order:** resolve the physical selections and materialize the
numeric decks; review P1's and P2's bounded cost and intended observation windows;
run the two precursors only after simulation authorization; assess the criteria
above; then finalize the production grid, timestep, duration, diagnostics, and
cost for a separate production decision. **Historical 2026-09-11 design-only boundary:** no simulation was started by
that update. The 2026-09-25 user request and execution record above now authorize
and record the bounded P1 segment; P2 and full production remain unsubmitted.

---

## CRITICAL PATH TO §B.3 PRODUCTION (state at 2026-08-31)

### Done and validated

| item | status |
|---|---|
| Stage-A no-impact gate | ✅ `λ = 0.2184495545645199`, enthalpy `4.975384548799866e-10`, max speed 0 — **bit-identical in all 16 digits on all SIX rounds** (20198502, 20213564, 20215987, 20216530, 20268028, 20273751), including after the `Velocity.c` / `msolve_cg.c` / `msolve_direct.c` changes |
| B.1.4 sediment clip | ✅ resolved (round 4, CH no-flux mask). `net = 0.000000e+00` in 2-D and 3-D |
| Meltwater tracer | ✅ **3.25e-3 in 3-D with a resolved grain**, bounded (was 9.3e-2) |
| Motion lock | ✅ exact — `U = 0`, `ΔY = 0` to machine zero, 2-D and 3-D |
| Interface-triggered release | ✅ fires in 2-D (`t = 19.03`, reproduces the recorded baseline exactly) and **in 3-D (`t = 311.38`)** |
| Impulse-free handoff | ✅ 0.24–0.52× `a_free` (bar < 2×) |
| `F_release` calibration | ✅ **0.709**, reconfirmed against the non-solid normalisation |
| Budgets | ✅ enthalpy 2.7e-4, salt 0, `C_S` drift exactly 0 |
| Settling pathology | ✅ root-caused (`darcy_tau` 3,906× too stiff) and fixed — `U_t` −6.1e-3 → **−0.56** |

### Current dependency gates — updated 2026-09-07

The [closure report](YANG_ECCO_B2_CLOSURE_REPORT.md) is the controlling numerical
record. Historical sections above preserve failed hypotheses and earlier jobs.

| # | Item | Current evidence / remaining action |
|---|---|---|
| 1 | Corrected settling refinement | **Pass for transients:** genuine ice-free controls, d/dx=24→32 changes speed by 0.935% at t=6 and ≤1.530% over t=1–6. No terminal plateau or Faxén-based confinement verdict is claimed. |
| 2 | Isolate sediment/ice coupling | **Pass:** tau=1e-3 and negligible-Darcy arms are bit-identical; conventional IBM agrees to roundoff through settling/collision. Three control jobs timed out after these observations. |
| 2b | Default-off no impact | **Pass:** final Stage-A job 20487628 reproduces every gate metric exactly. Previous convective legacy-path check 20274565 remains recorded separately. |
| 3 | Ice rigidity | **Pass in completed coupled run:** tau=1e-3 + hard threshold; worst sampled bulk/domain speed ratio 0.583%, within 1%. Historical 0.060% was one separate diagnostic, not a universal bound. |
| 3b | Phantom ice and periodic geometry | **Closed:** actual-ice transport preserves zero ice; conservative moving-mask remap; minimum periodic particle distance. Manufactured remap is mass-conservative and bit-identical on 8/16 ranks. |
| 4 | Clean coupled B.2 rerun | **Five checks pass, 20478206:** exact lock, smooth release, conserved tracer, settling/wall collision, heat/salt budgets. Version-2 restart and error guards also pass. |
| 5 | Select physical ECCO scenario / EOS | **Held for user approval and physical choices.** No real scenario or scenario EOS fit has been run. |
| 6 | Scenario-specific Cn, melt band, liquidus slope, mesh and cost | Follows the selected physical scenario; synthetic-debug values are not dimensional production defaults. |
| 7 | Timestep validation | **One test pending: 20487690**, d/dx=16, max_dt=0.005 with 20 ramp steps versus baseline 0.01/10. Coarse half-dt run changed physical release duration and is not a clean convergence test. |
| 8 | In-solver profiles and integrated dissipation | **Implemented and verified:** default-off ECCO_PROFILES, every 10 steps plus field output; liquid-weighted scalar fluxes, dissipation, phase/scalar budgets. Offline heat flux integration and source-aware interpretation; no unsupported mixing-efficiency claim. |

The five B.2 debugging checks are satisfied. Complete item 7 before declaring the
pre-scenario technical validation finished. Real ECCO simulation remains an
explicit user-approval boundary. Every closure run uses synthetic diagnostic
parameters, including passive nonzero salt for conservation checks.

### Standing methodological notes from this campaign

1. **Never key a locality argument on a tanh tail** — use the mesh (B.1.4 round 1).
2. **A silent fallback hides a broken fix** — report what a local step actually placed
   (B.1.4 round 1).
3. **Check the deck can express the failure you are testing for** — `StageB_release2d` and
   `StageB_release3d` both had no-flux thermal BCs and could never release.
4. **Agreement between runs is not validation against physics** — three runs reproduced a
   settling velocity that was 200× wrong, and matched the recorded baseline exactly.
   Validate against an independent expectation (a correlation, an exact solution), not
   against your own history.
5. **Measure the quantity your fix is premised on, before designing the fix.** The
   `darcy_tau` "window" at 1e-3 was derived from an assumed `φ_s ≈ 0.02` near the grain.
   The real value was 0.15–0.30 — and rose as damping weakened, which was the tell that it
   was a feedback. One cheap measurement first would have skipped a whole cycle.
6. **Test a fix cheaply before committing a long run** — the first 3-D submission burned
   12 h × 64 ranks against an unvalidated fix.
7. **A kernel validated at one liquidus slope is not validated at another** — the band's
   volume-averaged salinity costs ~0.5×`liquidus_slope` in interface temperature: harmless at 0.014
   (Yang, B.2), sign-flipping at ~1 (P1a). Re-gate the melt kernel whenever the regime changes.
9. **The Yang salt operator needs sediment treated as non-ice.** Its capacity 1/(F+δ) is
   about 10⁶ at C_L = 0 support cells beside water, and the moving grain made it diverge (P1b
   job 20973340). build_v3 uses C_L + C_S under `VOF_IBM`; it is under test in dev3 (21006068).
8. **Check that every division by a phase fraction is safe inside the resolved solid.**
   In B.3 P1b, liquid-referenced salinity s/(C_L+δ) was validated in fluid-only gates but
   reaches ~10³ inside a resolved grain. The symptom was a pressure loop stuck at a
   round-off-quantized divergence. A build that prints the divergence location found it
   in one short job.

---

## RESOURCE STATUS AND FORWARD REQUIREMENT (as of 2026-08-31)

**Allocation `phy250167`: 386,853 SU limit, 309,843 used, 77,010 SU remaining (20 %).**

### What the Stage-B campaign has cost so far

Whole ECCO campaign, 2026-08-28 → 2026-08-31: **3,425 SU — 1.1 % of all usage to date.**
The 309,843 SU already spent is overwhelmingly the Stage-A / Yang work, not this.

| category | SU | share of campaign |
|---|---:|---:|
| 3-D B.2 shakeout (5 runs incl. 1 timeout, 1 cancelled) | 3,188 | 93.1 % |
| 2-D B.1.4 A/B arms (13 runs) | 235 | 6.9 % |
| Stage-A no-impact gates (5 runs) | <1 | 0.0 % |
| settling gate (2 runs, in flight) | <1 | 0.0 % |
| **total** | **3,425** | |

Stated plainly: the four-round B.1.4 debugging was **cheap** (235 SU of 2-D arms); what
cost real time was the 3-D shakeout, and of its 3,188 SU roughly **768 SU (24 %) was
wasted** on the timed-out first run and on two runs whose thermal BC made the release
unreachable. Both causes are now fixed and recorded.

### Measured cost basis

**2,763 cell-steps per core-second**, from job 20218407 (655 k cells × 5 600 steps on 64
ranks in 5.76 h). This is a **pessimistic floor**: the shakeout runs a 22³ per-rank
subdomain and is communication-bound, so a production decomposition (32³–64³ per rank)
should recover ~3×. Quoted below as "measured" and "tuned".

| run | cells | steps | measured SU | tuned SU |
|---|---:|---:|---:|---:|
| B.3 baseline `d/Δx = 24`, `dt = 1.25e-2` (CFL 0.3) | 35.4 M | 8 000 | 28,463 | 9,488 |
| B.3 baseline, `dt = 2e-3` (**as written in the deck**) | 35.4 M | 50 000 | 177,894 | 59,298 |
| B.3 hero `d/Δx = 48` | 283 M | 16 000 | 455,408 | 151,803 |

> **The deck's `max_dt = 2e-3` costs a factor of 6 for nothing.** CFL 0.3 at `Δx = d/24`
> with `|u| ~ 1` allows `dt = 1.25e-2`, and the melt-source limit is `Pe_T·ε·4Cn = 0.034`,
> so CFL governs. The 3-D shakeout ran stably at `dt = 1e-2` throughout. **Confirming the
> true limit is the single highest-leverage thing to do before requesting time.**

### Forward requirement

| scenario | low (tuned) | high (measured) |
|---|---:|---:|
| Finish Stage B only (clean B.2 + settling gate + iterations) | 8,000 | 8,000 |
| + 1 baseline production run | 17,488 | 36,463 |
| + 3-run reduced sweep | 36,463 | 93,389 |
| **+ 7-run one-at-a-time sweep** | **74,414** | **207,241** |
| + 7-run sweep **and** 1 hero convergence run | 226,216 | 662,648 |

**Against 77,010 SU remaining:** enough to finish Stage B and do **one** baseline
production run comfortably, or a 3-run reduced sweep at the optimistic end. **Not** enough
for the 7-run sweep §B.3 specifies, and nowhere near enough for the hero run.

Note the full factorial in §B.3 (`Ga`×`St`×`Sc`×EOS = **36 runs**) is infeasible at any
of these numbers and should be replaced in the write-up by the one-at-a-time design costed
above.

### Recommended request

- **Minimum to deliver a publishable Stage-B result:** ~**110,000 SU** — finishes Stage B,
  the 7-run sweep at tuned throughput, and ~30 % contingency.
- **Recommended:** ~**300,000 SU** — the same plus the `d/Δx = 48` hero run, which is what
  makes the baseline defensible as a converged DNS rather than a single-resolution result.

Both assume the `dt` question above is settled first; if `max_dt = 2e-3` turns out to be
required, multiply the production lines by **6**.

### Caveat on the estimate

These extrapolate a **freshwater, `Sc = 7`** shakeout. Production runs `Sc = 7` **and
`Sc = 70`**, and the `Sc = 70` salinity field has an 8.4× thinner Batchelor scale on the
same single grid. That does not change the cell count (the grid is set by `d/Δx ≥ 24` for
the IBM) but it **does** raise the scalar-solver iteration count, and the estimate above
does not capture it. Treat the `Sc = 70` arms as the least certain line items.

---

## Quick reference — new flags, inputs, init types

| New flag | Stage | Purpose |
|---|---|---|
| `EOS_NONLINEAR` | A | Roquet quadratic buoyancy `b(θ,s)` |
| `CONC_VOF_PHASEWEIGHTED` | A | κ(F) for T, D_S·F for S (salt only in liquid) |
| `PHASE_CHANGE` | A | Stefan source in CH + latent sink in T + liquidus |
| `ICE_PENALIZATION` | A | Darcy damping `−(φ_s/τ)u` to keep ice rigid |
| `VOF_IBM` + `LAG_PARTICLE_RESOLVED` + `PARTICLE_RELEASE` | B | sediment particle, interface-triggered release |
| `VOF_DIFFUSE_SEDIMENT_CH_MASK` | B | **the B.1.4 fix.** No-flux (90°) condition at the grain surface in the CH biharmonic: identity rows for sediment cells + cancellation of stencil faces reaching into sediment. **Required in every ECCO deck with a resolved grain.** `#error` under `AXISYM_RZ` |
| `VOF_DIFFUSE_SHELL_FRACTION_NONSOLID` | B | release criterion normalised by non-solid shell volume; removes the `φ_liq` ceiling artifact. **Requires recalibrating `F_release`** |
| `VOF_DIFFUSE_SEDIMENT_CLIP_RESTORE` | B | ~~rim-local restore of clipped mass~~ **SUPERSEDED by `..._CH_MASK`, do not use** (§B.1.4-fix rounds 1–2) |
| `VOF_DIFFUSE_SEDIMENT_CLIP_AUDIT N` | B | log-only `SED_CLIP_AUDIT` lines every `N` steps: sediment trim, un-placed residual, and the fraction of the melt budget involved |

| New input block | Keys |
|---|---|
| `[eos]` | `eos_q, eos_betaT, eos_betaS, eos_Tmd0, eos_Tmd_slope` |
| `[phase_change]` | `stefan, T_melt, liquidus_slope, melt_band_eps` |
| `[vof]` | `darcy_tau`, `vof_slab_x0` (init-27 interface location; default 0.9) |
| `[conc]` | `kappa_ice_ratio` (per-field array), `theta_ice` (+ reuse `cbd0/cbd2/cbd5`, `conc_init_type`) |
| `[particle]` | `F_release` |

| New init type | Field | Meaning |
|---|---|---|
| VOF 27 | F | vertical ice slab (Yang `vof_slab_x0=0.9`; Stefan A.3 `0.05`; salty Stefan A.3b `0.5`) |
| VOF 29 | F | top ice layer + embedded sphere + water (ECCO) — **was 28; renumbered
2026-07-30 because `case 28` in `src/Cart3d.c:1150` is already the Favier melting-RB
horizontal solid layer. Using 28 for ECCO would have silently redefined a validated
testcase.** |
| Conc 30 / 31 | T / S | Yang: T step ; S linear-in-y, 0 in ice (A.3b: `cbd2 = cbd5` ⇒ uniform s) |
| Conc 32 / 33 / 34 | T / S / tracer | ECCO: two-layer T ; stratified S ; meltwater tracer = 0 |
| Conc 36 / 37 / 38 | T / S / tracer | ECCO B.3 P1b: erf diffusive sublayer below the VOF-28 ice (keys `sublayer_ell_T`, `sublayer_ell_S`, `sublayer_s_interface`; far values `cbd0`/`cbd5`, interface T `theta_ice`); tracer = meltwater fraction  Optional `sublayer_y0` anchors the profiles independently of the ice (no-ice twins). With `yang_salt_transport = 1`, init 37 writes the liquid salinity unmasked. |

**Two slopes that look alike but are not:** `eos_Tmd_slope` (EOS density-maximum shift, from `cS`)
≠ `liquidus_slope` (freezing-point depression, from the seawater freezing slope `m`). Keep them
distinct in code and inputs.
