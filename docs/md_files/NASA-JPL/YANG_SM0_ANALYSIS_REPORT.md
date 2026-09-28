# Yang `Sm = 0` freshwater-reference simulation analysis

> **VERDICT REVERSED 2026-07-27 — the "freshwater PASS" was scored against a bad reference.**
> This report's central result, `t½,0 = 60.2866` vs "Yang-inferred ≈ 59 (+2.2 %)", used a value
> *derived* from figure 4(a) on the assumption that figure 4(a) is the `ΔSv = 5` series. It is the
> **`ΔSv = 0`** series, and it contains the `Sm = 0` curve **directly**: `t½ = 107.42`. PARTIES is
> therefore **1.78× too fast** in freshwater, not within 2.2 %. The 59 was `0.5 × 118` from the
> broken chain, and PARTIES matching it was a coincidence of that construction.
>
> Everything else here stands: the run's numerics, conservation, morphology, the matched-transport
> audit, and the salty/fresh **ratio** result (`0.3352` vs `0.497`) are unaffected, because the
> ratio never depended on the misidentification. What changes is that the freshwater case is now
> the campaign's **largest** discrepancy rather than a passed gate — and, containing no salt, its
> cheapest. Provenance and corrected targets:
> [YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md); the follow-up ladder
> is `ECCO_TESTS/StageA_Yang/Yang_fresh_convergence/`.

**Decision: PASS the freshwater-baseline discriminator; FAIL the paired Yang salty/fresh
melt-rate ratio.** The completed `S \equiv 0` simulation gives
`t_{1/2,0} = 60.2866` and `\bar f_0 = 0.0165874`, within about 2.2 % of the
`t_{1/2,0} \approx 59` value inferred from Yang et al.'s published
`\bar f/\bar f_0 \approx 0.5` and the digitized `Sm = 5`, `\Delta S_v = 5`
half-melt time `t_{1/2,5} \approx 118`. The PARTIES pair instead gives
`\bar f_5/\bar f_0 = 0.3352`, only 67 % of the published normalized ratio.
The new run therefore validates the freshwater reference path and strengthens the conclusion
that the remaining Stage-A discrepancy is salt-specific double-diffusive heat delivery, or
uncertainty in the published/digitized Yang reference, rather than the basic Stefan/interface
machinery.

**Controlling follow-up 2026-07-21:** the paired `Sm=5` transport audit has been returned and
the precommitted rule has been applied. Salt recovery is resolved by tens of cells, leakage is
negligible, and there are no starved/refreezing rows; therefore **do not run `Sc=500`**. A
subsequent cheap Le=100 equation gate also found that mapping Yang's finite-interface salt
operator onto the existing PARTIES Cahn--Hilliard phase field worsens the analytic salinity
profile and front constant. Do not refine that gate or launch another PARTIES Yang production
case. See [YANG_SALT_OPERATOR_GATE_REPORT.md](YANG_SALT_OPERATOR_GATE_REPORT.md) for the
controlling result and next path.

Analysis date: 2026-07-21. Reference: Yang, Howland, Liu, Verzicco & Lohse,
*Ice melting in salty water: layering and non-monotonic dependence on the mean salinity*,
J. Fluid Mech. 969 (2023) R2
([DOI](https://doi.org/10.1017/jfm.2023.582),
[arXiv:2302.02357](https://arxiv.org/abs/2302.02357)). Read this report with the
[salty production analysis](YANG_PRODUCTION_ANALYSIS_REPORT.md) and the
[melting-RB campaign analysis](MELTING_RB_ANALYSIS_REPORT.md).

## 1. Run completion and provenance

| Item | Value |
|---|---|
| Directory | `/anvil/scratch/x-mjalabert/Yang_Sm0_07182026` |
| SLURM job | `19359337`, 2 nodes / 256 ranks, `wholenode` |
| State | `COMPLETED`, exit code `0:0` |
| Runtime | 39 h 05 min 28 s; 2026-07-19 05:25 to 2026-07-20 20:30 EDT |
| Grid | `1440 x 1440 x 1`, `TWOD_CARTESIAN`, periodic storage direction `z` |
| Evolution | 400,019 iterations; `t = 200.000006`; `dt` from `1e-4` to the `5e-4` cap |
| Outputs | 101 complete `Data_*` and 101 complete `Resume_*` files, indices 0--100 |
| Executable SHA-256 | `aac878bfb6cdfe8af198685f80144388c25838231befb5ffe4865c03ddf65d38` |
| Input SHA-256 | `2992ade15971aeb3ba0df9d96d152c262966ef9fe15ade96e56fc2c1429f87ac` |
| Job-script SHA-256 | `50c665e2aecb8a413890282bdc501ec5dbccce4a83cf6e99712b5276ba2c54e9` |

The archived executable is byte-identical to the repository executable at analysis time.
Relative to the `Sm = 5`, `\Delta S_v = 5` production input, the physical input changes are
only `cbd2 = cbd5 = 0`, so the salinity field is identically zero. The retained salt-dependent
EOS and liquidus coefficients multiply `s = 0` and therefore contribute exactly zero. This is
the correct no-salinity normalization run `\bar f_0`.

All HDF5 calculations below use only the physical cells `[:nz, :ny, :nx]`; the duplicated
high-side PARTIES storage planes are excluded. The measured initial ice volume is
`V_0 = 0.10000000000000003`.

## 2. Reference interpretation

Yang et al. define `\bar f = 1/t_{1/2}` and normalize it by the melt rate without salinity,
`\bar f_0`. Their final-paper figure 4 is figure 3 in the arXiv letter. Two distinct quantities
must not be conflated:

1. Figure 4(a) is stated to hold `\Delta S_v = 5 g kg^-1`; its curve labelled `Sm = 0`
   is therefore not safely interchangeable with the `S \equiv 0` normalization used for
   `\bar f_0` in figure 4(b). The paper does not publish the raw `S \equiv 0` time series.
2. For the reference `Sm = 5`, `\Delta S_v = 5` point, the project digitization gives
   `t_{1/2,5} \approx 118`, while figure 4(b) gives `\bar f_5/\bar f_0 \approx 0.5`.
   These imply

   `t_{1/2,0} = (\bar f_5/\bar f_0) t_{1/2,5} \approx 0.5 x 118 = 59`.

Because both source values are read from plotted curves, `59` is an inferred target with
figure-reading uncertainty, not author-supplied source data. The outstanding data-request
draft in [LOHSE_DATA_REQUEST_EMAIL.md](LOHSE_DATA_REQUEST_EMAIL.md) explicitly asks for the
raw freshwater normalization and figure-4 data.

## 3. Melt-rate results

| Quantity | PARTIES | Yang target/inference | Status |
|---|---:|---:|---|
| Freshwater half-melt time `t_{1/2,0}` | **60.2866** | approximately 59 | **PASS** (+2.18 %) |
| Freshwater rate `\bar f_0` | **0.0165874** | approximately 0.01695 | **PASS** (-2.14 %) |
| `V(200)/V_0`, freshwater | **0.197875** | no raw `f_0` curve published | diagnostic only |
| Salty production half-melt time `t_{1/2,5}` | 179.8484 | approximately 118 | FAIL, as previously reported |
| Paired `\bar f_5/\bar f_0 = t_{1/2,0}/t_{1/2,5}` | **0.335208** | approximately 0.5 | **FAIL** (0.670 of target) |

Selected freshwater history:

| `t/t_ff` | `V/V_0` | max speed | relative enthalpy drift |
|---:|---:|---:|---:|
| 40 | 0.630589 | 0.1270 | -8.00e-6 |
| 60 | 0.501679 | 0.1221 | -9.32e-6 |
| 80 | 0.397293 | 0.1015 | -5.28e-5 |
| 100 | 0.326738 | 0.0609 | -1.91e-4 |
| 120 | 0.281322 | 0.0404 | -2.72e-4 |
| 160 | 0.227481 | 0.0229 | -3.24e-4 |
| 200 | 0.197875 | 0.0134 | -3.55e-4 |

The half-melt discriminator is complete by `t = 60.3`; extending to `t = 200` establishes
the late morphology and exposes a separate wall-melt-through budget issue, but does not change
the rate-ratio decision.

## 4. Numerical health and conservation

| Check | Result | Status |
|---|---:|---|
| Log completion | 256 coordinated `Exit here.` lines; no failure token or NaN/Inf | PASS |
| Projection | maximum divergence `1.1e-9` | PASS |
| Solver iterations | CH <= 17; velocity <= 9; concentration <= 4; HYPRE <= 9 | PASS |
| Phase bounds | `F_min = 0`, `F_max = 1` over all outputs | PASS |
| Out-of-bounds counters | zero in both scalar files over the entire run | PASS |
| Salinity | exactly zero in every analyzed cell/output | PASS |
| Maximum fluid speed | 0.136314 | active convection |
| Maximum speed in `F < 0.1` | 1.686e-4 | adequately penalized |
| Pre-half two-output enthalpy change | maximum `5.996e-7` through `t = 60` | PASS |
| Pre-half cumulative enthalpy drift | `-9.325e-6` at `t = 60` | negligible for `t_{1/2}` |
| Full-run cumulative enthalpy drift | `-3.545e-4` at `t = 200` | late warning |
| Largest two-output enthalpy change | `1.434e-5` at `t = 84` | late warning |

At half melt, the cumulative enthalpy drift is equivalent to only `0.0105 %` of the initial
ice volume. At `t = 200`, it is equivalent to `0.399 %` of the initial ice volume. It is much
too small to create the melt-rate result, but it exceeds the desired late-time budget criterion.

The timing localizes the issue. The first output with fully melted rows appears at `t = 72`;
the largest two-output budget change follows at `t = 84`. Before local wall melt-through,
the two-output budget closes below `1e-6`. The late drift is therefore most consistent with
phase-field/interface annihilation where the retreating ice meets the east wall, rather than
with a bulk latent-heat conversion error. It should be audited independently, but it does not
invalidate `t_{1/2,0}`.

## 5. Freshwater morphology

There is no salt stratification and therefore no Yang/Huppert--Turner layer-spacing target for
this run. The morphology instead shows the expected thermally driven asymmetry and progressive
loss of the upper ice:

| Time | bottom 2 % ice thickness | mid-height thickness | top 2 % thickness | fully liquid rows |
|---:|---:|---:|---:|---:|
| 40 | 0.07963 | 0.06593 | 0.03819 | 0 |
| 60 | 0.07539 | 0.05461 | 0.01304 | 0 |
| 72 | 0.07378 | 0.04858 | approximately 0 | 65 |
| 100 | 0.07191 | 0.03592 | 0 | 388 |
| 160 | 0.06856 | 0.00084 | 0 | 730 |
| 200 | 0.06642 | 0 | 0 | 847 (58.8 %) |

The remaining ice is a bottom-attached wedge by `t = 200`, while the circulation decays as the
thermally active interface area shrinks. This is qualitatively reasonable for the nonlinear
freshwater EOS and contains no layer-forming salt mechanism.

## 6. Interpretation for Stage A

The result updates the causal picture as follows:

1. **Freshwater Yang-geometry melting passes.** The same geometry, nonlinear thermal EOS,
   phase-change kernel, penalization and `Ra_T = 1e7`, `Pr = 10` parameters reproduce the
   inferred freshwater half-melt time within a few percent.
2. **The salt-free validation is now direct, not only transferred from melting-RB.** Together
   with the Favier St = 0.1 result and flat-wall RB Nusselt check, this rules out a systematic
   salt-free convective-melt deficit in PARTIES.
3. **The paired Yang response remains wrong.** PARTIES gives `\bar f_5/\bar f_0 = 0.335`,
   versus approximately `0.5`; the salty case is suppressed about one-third too strongly in
   normalized-rate terms.
4. **Do not tune the Stefan source or Darcy penalty from this result.** The quiescent gates,
   in-situ interface audit, Darcy and CH-mobility null branches, melting-RB budgets, and now
   the direct freshwater half-melt time all point away from those components.
5. **Stage A remains open.** The morphology criterion passes in the salty case, and the
   freshwater baseline passes here, but the required salty/fresh response has not been matched.

## 7. Next steps

### 7.1 Resolve the reference before another production run

Send the prepared [data request](LOHSE_DATA_REQUEST_EMAIL.md) and obtain, if possible:

- the raw `S \equiv 0` `\bar f_0` or `V(t)` series;
- the source data behind final-paper figure 4(a) and 4(b);
- confirmation of how the `Sm = 0`, `\Delta S_v = 5` label/profile was constructed;
- the sign convention for `\Delta S_v`.

The present inference is internally consistent and the PARTIES freshwater match is strong, but
author data would remove the largest remaining reference uncertainty before expensive reruns.

### 7.2 Matched transport audit -- paired audit complete 2026-07-21

The allocation-free audit is implemented in
[`analyze_yang_matched_transport.py`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/analyze_yang_matched_transport.py).
It scans only `time` and `F` in all outputs, reads full fields only at the four selected states,
and reconstructs the exact instantaneous PARTIES Stefan source

`m = St/(Pe_T eps) (theta - theta_L) F(1-F)/(sqrt(2) Cn)`.

It also records a liquid-side thermal/salinity 90 % recovery-length proxy, phase-weighted salt
diffusive flux at `F = 0.5`, interface velocity and buoyancy, salt leakage, starved/refreezing
row fractions, and twelve vertical melt-delivery bins. The thermal recovery length is measured
relative to the value 128 cells into the liquid; it is a repeatable paired diagnostic, not an
assertion that the convecting temperature field has a classical monotone boundary layer.

Freshwater results are:

| Target `V/V0` | Selected `V/V0` | `t` | reconstructed `<m>` | `<m>` / finite-difference rate | interface superheat | median `delta_T,90/dx` | interface `<|v|>` | row-melt CV | starved rows |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.8 | 0.808988 | 18.0 | 9.4982e-4 | 1.0060 | 0.01978 | 106.98 | 2.914e-3 | 0.347 | 7.64 % |
| 0.7 | 0.705410 | 30.0 | 7.9391e-4 | 1.0027 | 0.01664 | 106.98 | 2.199e-3 | 0.482 | 14.10 % |
| 0.6 | 0.602858 | 44.0 | 6.8118e-4 | 1.0024 | 0.01434 | 106.17 | 1.617e-3 | 0.579 | 17.71 % |
| 0.5 | 0.501679 | 60.0 | 5.9109e-4 | 1.0012 | 0.01248 | 98.07 | 1.147e-3 | 0.636 | 17.64 % |

The 0.12--0.60 % source/rate closure independently validates the source reconstruction and
the saved-field slicing. The freshwater baseline shows the expected decline in interfacial
speed and mean melt delivery, with increasing vertical heterogeneity as the bottom-attached
ice shape develops. Salt recovery and salt-flux values are correctly absent/zero for `S = 0`.

The salty raw archive, `/bigscratch/mjalabert314/Yang_production`, was processed data-locally and
only the small products were returned. The executed command was:

```bash
python3 "$HOME/PARTIES/PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/analyze_yang_matched_transport.py" \
  --run sm5=/bigscratch/mjalabert314/Yang_production \
  --baseline-matched-csv \
    "$HOME/PARTIES/PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/matched_transport_audit_sm0/yang_matched_transport_matched.csv" \
  --output-dir /bigscratch/mjalabert314/Yang_production/matched_transport_audit
```

The returned audit verified all 101 outputs through `t=200.00000638`, the required schema,
matched snapshots, and physical-domain slicing. Its controlling findings are:

- reconstructed source agrees with finite-difference melt rate within 0.53 %;
- `Sm5/Sm0` melt-rate ratios are 0.25--0.37, heat-flux ratios are 0.40--0.47, and superheat
  ratios are 0.25--0.36;
- every salinity profile is usable, no profile has `delta_S,90 < 4 dx`, and median recovery
  lengths are 27.7--98.4 cells;
- salt overlap/leakage indices are only `1.25e-4`--`1.92e-4`; and
- there are no starved or refreezing rows under the specified definitions.

The low melt is therefore consistent with insufficient convective heat delivery, not a
numerically unresolved salt film or anomalous salt leakage.

### 7.3 Allocation-minimizing decision after the salty audit -- applied

**Outcome: condition 3 applies. Do not run `Sc=500`.** The following rule is retained as the
precommitted provenance for that decision:

Do not start a new simulation until the salty CSV exists. Compare salty/fresh values at each
requested ice fraction, emphasizing reconstructed `<m>`, liquid-side heat-flux proxy,
interface vertical speed, row-melt starvation, vertical melt-delivery bins, `delta_S,90/dx`,
salt-flux direction, and the solid-salt leakage index.

Use the following precommitted decision rule:

1. If the salty recovery layer is marginal (`delta_S,90 < 4 dx` over at least 25 % of usable
   interface rows), or if anomalous phase-weighted salt flux/leakage is vertically correlated
   with the heat-starved rows, run exactly one `Pe_S = 5e5` (`Sc = 500`) branch from the existing
   salty `t = 20` checkpoint through `t = 60`. Keep the 1440 grid and every other input fixed.
2. At `t = 60`, the existing volume gap is 0.07763 in `V/V0`. Continue this direction only if
   the branch removes at least 20--25 % of that gap (a decrease of 0.0155--0.0194 in `V/V0`),
   while salt drift remains below `1e-8` and the layer morphology remains coherent. A response
   below 5 % of the gap (less than 0.0039 in `V/V0`) is a stop result; do not add another `Sc`.
3. If the salt film is resolved, leakage is negligible, and its transport metrics do not
   explain the missing heat delivery, do not run the `Sc = 500` branch. Treat the remaining
   discrepancy as reference/method-form uncertainty and use the existing data, paper, and AFiD
   implementation for further reconciliation.
4. Do not repeat Darcy or CH-mobility branches, perform a parameter sweep, run the separate
   late-wall cleanup case, or launch another `t = 200` production simulation while allocation
   is constrained. Only a successful short discriminator can justify one final full run.

### 7.4 Audit late wall melt-through separately

Use a reduced fresh-water case that reaches the east wall cheaply and close
`d<theta + F/St>/dt` against all boundary fluxes and the integrated latent source across the
first local melt-through. Inspect conservative clipping/mass restoration and the last diffuse
interface cells at the no-flux wall. This is a cleanup/robustness task; it is not the cause of
the Yang rate-ratio miss.

### 7.5 Closure rule — REVISED 2026-07-27

**Superseded text (2026-07-21):** *"An accurate Yang reproduction should use the complete
AFiD-MuRPhFi Allen--Cahn/Hester phase--salt method, or port and validate that complete coupled
update in PARTIES before more 2-D allocation is spent."* That route is **withdrawn as off-goal** —
the benchmark must be reproduced with PARTIES, and no evidence identifies the conservative
Cahn--Hilliard formulation as the defect (the Le=100 gate in fact ranks it **above** the mapped
Yang salt operator against an exact similarity solution).

**In force instead:** do not launch another `1440²` `t → 200` production run yet, but do proceed
with the PARTIES-only ladder in
[YANG_ECCO_IMPLEMENTATION_ROADMAP.md](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) §A.5 — free reference
audit and clip diagnostics (§A.5.1), the missing `720²/1440²` grid-and-interface-width convergence
pair at `Sm=5` (§A.5.2, ≈ 375 core-h), the `f̄(Sm)` curve-shape test at `Sm = 10, 15` (§A.5.3), and
interface-model variants screened in the 1-D Le=100 sandbox with an Allen--Cahn port ranked last
(§A.5.4). The quantitative gates are unchanged: a full run is authorized only after ≥ 25 % of the
`t=60` volume gap is closed (`V(60)/V₀ ≤ 0.7549`) with healthy conservation and morphology, and
Stage A closes only when the salty/fresh rate ratio, absolute salty volume history, `0.287H` layer
scale, and all conservation/solver gates pass together. Until then, Stage B remains blocked.

## 8. Saved artifacts

- Full 101-output physical-domain time series:
  [`yang_sm0_final_timeseries.csv`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/yang_sm0_final_timeseries.csv)
- Machine-readable final metrics:
  [`yang_sm0_final_metrics.csv`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/yang_sm0_final_metrics.csv)
- Summary figure:
  [`yang_sm0_final_analysis.png`](../../figures/yang_sm0_final_analysis.png)
- Matched-state transport audit:
  [`yang_matched_transport_matched.csv`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/matched_transport_audit_sm0/yang_matched_transport_matched.csv)
- Full lightweight volume scan and twelve-bin vertical diagnostics:
  [`yang_matched_transport_scan.csv`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/matched_transport_audit_sm0/yang_matched_transport_scan.csv) and
  [`yang_matched_transport_vertical_bins.csv`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/matched_transport_audit_sm0/yang_matched_transport_vertical_bins.csv)
- Reusable matched-audit implementation:
  [`analyze_yang_matched_transport.py`](../../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/analyze_yang_matched_transport.py)
- Raw run: `/anvil/scratch/x-mjalabert/Yang_Sm0_07182026`
