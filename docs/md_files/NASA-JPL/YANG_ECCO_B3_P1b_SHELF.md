# B.3 P1b — released fine-sand grain beneath an Antarctic ice-shelf base

**Status (2026-09-28):** design fixed, code change built, pre-flight job
submitted; see [§8 Run record](#8-run-record). This file is the controlling
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
- `liquid_referenced_salinity = 1`, so the liquidus and EOS read s/(F+δ) and the
  band starts in equilibrium;
- `yang_salt_transport = 0`;
- this choice is conditional on the 1-D gate below.

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
| Salt drift | ≤ 10⁻⁸ relative |
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

*Filled in as jobs complete.*

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
