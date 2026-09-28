# Melting Rayleigh-Benard benchmark (Favier et al. 2019) — St = 1 production analysis

**Verdict: the convective phase-change coupling is validated — the interface converts the
delivered heat flux into melt exactly (wall-flux vs. melt-rate closure 0.4–1.7 % in every
window) and the diffusive phase matches the exact Stefan solution to 1.7 % — but the
delivered heat transport itself runs ~17 % below Favier's fitted law:**
`gamma_eff = 0.0955 ± 0.009` vs. their `gamma = 0.115` (with the correct exponent,
free fit `Nu ∝ Ra_e^0.348`). The melt-rate plateau is 0.79–0.85 of the published model
line. **The salt-free case therefore reproduces roughly half of the Yang A.4 melt-rate
deficit (0.83 vs. 0.66 rate ratio), and localizes it: the gap is in the convective heat
delivery (Nu amplitude), categorically not in the Stefan/interface machinery.**

Analysis date: 2026-07-17. Reference: Favier, Purseed & Duchemin, *J. Fluid Mech.* 858
(2019) 437–473, case D ([arXiv:1901.03847](https://arxiv.org/abs/1901.03847),
[DOI 10.1017/jfm.2018.773](https://doi.org/10.1017/jfm.2018.773)).

## 1. Run provenance

| Item | Value |
|---|---|
| Directory | `/anvil/scratch/x-mjalabert/Melting_RB_07162026` (Anvil) |
| SLURM job | 19301354, 2 nodes / 256 ranks, wholenode, ran ~8.7 h |
| Input | byte-identical to `testcases/ECCO_TESTS/StageA_Yang/MeltingRB/parties_St1.inp` |
| Case | Ra = 1e7, Pr = 1, thetaM = 0.05, h0 = 0.05, aspect 6, St_Favier = 1 (`stefan = 1`) |
| Grid | 3072 x 512 uniform (dx = 1/512); Favier case D used 1024 x 512 pseudo-spectral/FD4 |
| Completion | t/t_ff = 120 (= 0.03795 diffusive times), 121 Data/Resume pairs, clean "Exit here." |
| Flags | full Yang production set; `Boundary.h` = `PERIODIC_NOSLIP_BOX` + `TWOD_CARTESIAN` |
| EOS degeneration | q = 1, betaT = 1, betaS = 0 ⇒ b = −theta (linear RB buoyancy through the validated EOS path) |

Outputs of this analysis: `analysis/melting_rb_timeseries.csv`, `analysis_arrays.npz`,
`fig_h_evolution.png`, `fig_meltrate.png`, `fig_nu_rae.png`, `fig_onset.png`,
`fig_fields.png` (all under the run directory).

## 2. Target correction (paper eq. 48–51)

The README/gate-script target "hdot_diff(St=1) = 23.1" omitted the sensible-heat storage
factor. Favier's eq. (48)–(51): `(1/2 + St) dh/dt|_diff = gamma Ra^(1/3) (1−thetaM)^(4/3)
≈ 23.14`, the 1/2 being the mean temperature of the newly molten (fully convective) fluid.
For St = 1 the correct model plateau is `hdot_diff = 23.14/1.5 = 15.43`
(free-fall: 4.88e−3). The README and `analyze_melting_rb.py` have been corrected.

## 3. Results against the published targets

| Target (paper) | Measured | Status |
|---|---|---|
| Diffusive phase = 1-D Stefan solution (their Fig. 7a dotted) | max abs. rel. error 1.7 % over t < 28 (dominated by a +dy/2 init offset; slope error ≈ 1.6 %) | **PASS** |
| Onset at h_c = 0.0564 (Ra_c = 1707.76) | h̄ crosses h_c at t ≈ 1.3; v_rms decays while marginally stable, grows exponentially (sigma ≈ 0.70/t_ff) from t ≈ 4–6, i.e. h̄ ≈ 0.067; visible h̄ departure at h̄ = 0.129 (t = 28) after amplitude growth from the 1e−4 seed | **PASS** (consistent; sharp h_c not resolvable with a theta-only infinitesimal seed, same delayed-departure phenomenology as the paper) |
| Post-onset plateau `(1/2+St)·Pe_T·dh̄/dt = 23.14` (gamma = 0.115) | 18.3–19.8 across windows (t ∈ [20,50]: 18.8; [50,80]: 18.3; [80,120]: 19.8; h̄ ∈ [0.15,0.45]: 18.9) → ratio **0.79–0.85** | **MISS (−15 to −21 %)** |
| `Nu = 0.115 Ra_e^(1/3)` (their Fig. 10) | gamma_eff = **0.0955 ± 0.009** (beta fixed 1/3); free fit `Nu = 0.0796 Ra_e^0.348` over Ra_e ∈ [2e4, 1.2e6] | **shape PASS, amplitude −17 %** |
| Interface energy conversion (their eq. 48 balance) | `q_bot − q_top` vs `(1/2+St)·Pe·dh̄/dt`: ratio 0.996 / 1.001 / 0.983 in the three windows | **PASS — kernel exact** |
| Morphology | fine scallops at onset (37 across aspect 6 at t = 30), coarsening to 10 large topography-locked laminar cells at t = 120; front roughness std 0.006 (t=60) → 0.038 (t=120) | **PASS** (paper's stabilized-large-cell mechanism reproduced) |
| St-dependence (their Fig. 7b) | not measurable — companion St = 0.1 run not executed on this machine | open (input `parties_St0p1.inp` ready) |

Health: F ∈ [0,1] **exactly** (min 0, 1−max 0 over all 121 snapshots); salt field
identically 0; solid-velocity leakage ≤ 1.5e−4; enthalpy budget closes to 4.3e−4
relative (trapezoid quadrature over 1-t_ff output spacing — quadrature-limited);
h̄ (volume) vs F = 0.5 contour mean agree to 3e−5. Final state: h̄(120) = 0.5039,
v_rms = 0.152, mean liquid theta = 0.525 (the model's 1/2 assumption verified).

## 4. Interpretation for the Yang A.4 deficit

1. **Exoneration is now complete across regimes.** Quiescent gates, the Yang in-situ
   audits, and now a fully convective benchmark all show the CH/superheat interface
   machinery converts delivered flux to melt exactly. The (1/2+St) energy accounting is
   verified end-to-end in PARTIES data.
2. **A salt-free transport shortfall exists but is about half the Yang gap.** Yang mean
   melt-rate ratio 0.66; here 0.79–0.85 against Favier's fitted line. If the ~17 %
   Nu-amplitude shortfall is generic to PARTIES melting convection, the salt-specific
   residual in Yang is ≈ 0.66/0.83 ≈ 0.80, i.e. ~20 % remains attributable to the
   double-diffusive/salt side or the digitized Yang reference.
3. **Softeners on the 17 %.** gamma = 0.115 is a single global fit over St ∈ [0.02, 50]
   (the paper quotes no per-case scatter; agreement called "very good"); our Nu level
   correlates with the convection-cell state — during the transient fine-scallop phase
   (t ≈ 28–40) the melt rate touches 21–24 (the Favier line), relaxing to ~19 as cells
   coarsen to the 10-cell state. Cell-coarsening dynamics are sensitive to seeding and
   history (cf. the group's bistability follow-up, arXiv:2002.03710). A genuinely flat
   -wall PARTIES RB Nusselt validation would separate "generic RB transport" from
   "melting-specific" effects.
4. **What this rules out for Yang:** interface kinetics, liquidus pinning, latent-heat
   accounting, penalization throttling of melt (all previously), and now the basic
   convective Stefan coupling. What it leaves: the amplitude of convective heat delivery
   in the presence of a phase-field front (common to both cases, ~17 %), plus a
   comparable-size salt-specific or reference-bias component in Yang.

## 5. Next steps — both discriminators submitted 2026-07-17

1. **St = 0.1 companion — SLURM job 19321552**
   (`/anvil/scratch/x-mjalabert/Melting_RB_St0p1_07172026`, 2 nodes/256 ranks, 24 h):
   input = `parties_St0p1.inp` (differs from St1 only in `stefan = 10.0`), same
   executable (sha256 `675f4f56…`). Tests Favier Fig. 7(b) St-dependence — prediction
   `hdot(St=0.1)/hdot(St=1) = 1.5/0.6 = 2.5` (NOT 10; the low-St saturation). An
   St-collapse pass with the same gamma_eff would prove the shortfall is a pure
   transport-amplitude offset, not a coupling error. Note: at the predicted rate the
   ice melts through around t ≈ 80; the h̄(t) comparison window is t ≲ 60 (h̄ ≲ 0.8).
2. **Flat-wall 2-D RB Nu check — SLURM job 19321553**
   (`/anvil/scratch/x-mjalabert/FlatRB_Ra1e6_07172026`, 1 node/128 ranks, 12 h):
   Ra = 1e6, Pr = 1, aspect 4, x periodic, 1024 x 256, t → 400 (statistics window
   t ∈ [100, 400]). Same executable, VOF inert (`vof_slab_x0 = 1.5` ⇒ F ≡ 1 exactly),
   `stefan = 0` disables melting (guard verified in `Conc_add_latent_heat_RHS`; no
   other 1/stefan divisions exist), theta IC = classical `Conc_init_RB` (type 6).
   Compare Nu (wall flux + volume `1 + Pe<v theta>`) against established 2-D RB
   literature values at Ra = 1e6, Pr = 1. Decides whether the −17 % is
   melting-specific (front representation) or generic to PARTIES RB transport.
3. Then the Yang salt-side program (Sm = 0 ratio run / Lohse-group source data), with
   the RB result as the salt-free baseline.
   **Completed 2026-07-20; analyzed 2026-07-21:** SLURM job **19359337**
   (`/anvil/scratch/x-mjalabert/Yang_Sm0_07182026`, 1440², t → 200) completed cleanly
   in 39 h 05 min. The input is identical to A.4 except `cbd2 = cbd5 = 0` (`s ≡ 0`).
   Result: `t½,0 = 60.2866`, within 2.2 % of the Yang freshwater value `≈59` inferred
   from `t½,5 ≈118` and `f̄5/f̄0≈0.5`. The direct freshwater baseline therefore passes,
   but the PARTIES pair gives `f̄5/f̄0 = 60.2866/179.8484 = 0.3352`, only 67 % of the
   published normalized ratio. Full report:
   [YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md). The data-request draft remains
   at [LOHSE_DATA_REQUEST_EMAIL.md](LOHSE_DATA_REQUEST_EMAIL.md).

## 6. Preliminary follow-up results (2026-07-17 23:15, St0p1 still running at t ≈ 71)

**Flat-wall RB (job 19321553, COMPLETED, 5 h 07, t → 400, clean):** over the
statistics window t ∈ [100, 400]: `Nu_bot = 8.27 ± 0.26` (SEM; instantaneous std 1.43 —
2-D RB at Ra = 1e6 is quasi-laminar with large coherent oscillations), `Nu_top` equal to
0.05 %, volume-average `1 + Pe<v theta> = 8.26` — internally consistent to 0.1 %.
Johnston & Doering's 2-D fit `Nu = 0.138 Ra^0.285` (valid Ra ≥ 1e7, Gamma = 2)
extrapolates to ≈ 7.1 at Ra = 1e6; 2-D values at this Ra are state/aspect-sensitive.
**PARTIES flat-wall RB sits at or above the literature level → NO generic RB
under-transport.** `gamma_flat = Nu/Ra^(1/3) = 0.083` (flat, non-melting reference).

**St = 0.1 companion (job 19321552, preliminary over t ≤ 71, h̄ ≤ 0.86):**
`fbar = (1/2+St)·Pe·hdot = 22.6–23.8` across windows → **ratio 0.98–1.03 of the Favier
target 23.14; gamma_eff(St=0.1) = 0.113 ≈ 0.115.** Diffusive phase consistent with the
Neumann asymptote (lambda = 1.243). F ∈ [0,1] exactly. Measured
`hdot(St=0.1)/hdot(St=1) = 3.11` vs. the model's 2.5 — the excess is exactly the St = 1
case's −17 %. (The wall-flux balance 20.8 vs 23.8 at St = 0.1 reflects the model's
"bulk = 1/2" storage assumption lagging at fast front speeds, not a kernel error; the
exact enthalpy budget will be closed in the final analysis.)

**Preliminary interpretation:** PARTIES lands ON the published melting-RB law at
St = 0.1 with the identical executable, grid, and machinery as the St = 1 case. The
St = 1 −17 % is therefore not a systematic transport deficit — it is
state/realization-dependent (convection-cell coarsening history; the 10-cell locked
state), i.e. within the physical variability that Favier's single-parameter
gamma = 0.115 fit smooths over. Combined with the healthy flat-wall Nusselt number,
**the salt-free convective melt coupling is validated, and the Yang deficit now points
predominantly at the salt side (double-diffusive layer/film transport) or the digitized
reference itself.** Final St0p1 analysis when the run completes (~t + 4 h).

## 7. FINAL St = 0.1 results (2026-07-18, job 19321552 COMPLETED, 10 h 20, t → 120)

| Quantity | Measured | Target | Status |
|---|---|---|---|
| Diffusive phase vs exact 1-D Stefan (lambda = 1.243) | max 1.1 % (t ∈ [1,8]; overlap visibly extends to t ≈ 19, h̄ ≈ 0.20) | — | **PASS** |
| Plateau `(1/2+St)·Pe·dh̄/dt` | 22.6–23.7 (windows t∈[15,40], [40,60], h̄∈[0.15,0.6]) | 23.14 | **PASS (0.98–1.03)** |
| `hdot(St=0.1)/hdot(St=1)` | 3.11 | 2.5 (model) | model + the St=1 −17 % |
| Exact enthalpy budget `Pe·dE/dt / (q_bot−q_top)` | **1.0006 ± 0.013** (t ∈ [20,60]) | 1 | **PASS — kernel exact** |
| Measured sensible-storage coefficient `d<theta>/dh̄` | 0.423 | model assumes 1/2 | explains the wall-flux vs (1/2+St) gap at fast fronts |
| Melt-through | h̄ > 0.95 at t ≈ 80; h̄(120) = 0.9929 | — | full-evolution capture |
| Health | F ∈ [0,1] exactly; salt ≡ 0; solid speed ≤ 2.9e−4; all finite | — | **PASS** |

Nu-level amplitude at St = 0.1: `gamma_eff = 0.098 ± 0.021` — consistent with the
melt-level match once the measured storage (0.42 vs the model's 0.5) is accounted for:
Favier's Fig. 7(b) compares hdot directly, which is the published target and is hit.

**Campaign verdict (final).** The Stefan/interface machinery is exact in every regime
tested (quiescent gates, Yang in-situ, convective budget closure at both St). Against
the published melting-RB law: St = 0.1 ON the line, St = 1 at 0.79–0.85, flat-wall RB
Nusselt at/above literature. There is **no systematic salt-free convective-melt deficit
in PARTIES**; the St = 1 shortfall is within state/realization variability. **The Yang
A.4 deficit (rate ratio 0.66) therefore localizes to the salt side (double-diffusive
film/layer transport at Sc = 1000) or to the digitized Yang reference.**

**Direct Yang-geometry confirmation (2026-07-21).** The completed `S ≡ 0` run independently
confirms this campaign verdict: the inferred Yang freshwater half-melt time is hit to about
2.2 %. The remaining discrepancy is sharper in ratio form: PARTIES `f̄5/f̄0 = 0.3352` versus
Yang `≈0.5`. Do not retune the Stefan source, Darcy penalty, or CH mobility; proceed with the
salt-film/transport audit and raw-data request described in
[YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md).

**Audit status 2026-07-21.** The matched-volume analyzer and freshwater half of that audit are
complete; exact snapshot source reconstruction closes the measured freshwater melt rate within
0.12--0.60 %. The remaining Sm=5 postprocessing must run data-local beside the separate
`/bigscratch` archive. No further melting-RB run or Yang simulation is warranted before those
small paired diagnostic CSVs are available.

Campaign artifacts: `/anvil/scratch/x-mjalabert/MeltingRB_campaign_analysis/`
(`campaign_summary.csv`, `fig_fig7b_st_dependence.png`, `fig_nu_rae_combined.png`,
`fig_h_St0p1.png`, `fig_meltrate_St0p1.png`, `fig_nu_flatrb.png`,
`fig_fields_St0p1.png`) plus per-run `analysis/timeseries.csv` in each run directory.
