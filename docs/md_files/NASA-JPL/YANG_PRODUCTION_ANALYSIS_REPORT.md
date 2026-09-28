# Yang production simulation analysis

> **RE-SCORED 2026-07-27 — the targets in this report are wrong; the FAIL stands but reverses
> direction.** `t½ = 118` and `V(200)/V₀ = 0.32` were digitized from Yang figure 4(a), which is
> the **`ΔSv = 0`** series, not this run's `ΔSv = 5`. The correct target for `Sm=5, ΔSv=5` is
> `f̄/f̄₀ = 0.497` ⇒ `t½ ≈ 216`. PARTIES' `179.848` is therefore **1.20× too fast**, not "52.4 %
> later than the target". Every measurement, budget, morphology result and null branch below is
> unaffected. See [YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md).

**Decision: FAIL the Stage-A Yang benchmark, despite a qualified morphology pass and an
otherwise healthy run.** The required acceptance is conjunctive: the melt-front layer spacing
must agree with Yang/Huppert--Turner **and** the normalized ice-volume history must track Yang's
2-D reference curve. The layer spacing and top-thickness signature are reproduced, but
`t_1/2 = 179.8484` is 52.4 % later than the target `118`, and `V(200)/V0 = 0.45616` is 42.5 %
above the target `0.32`. Stage A therefore remains open.

Analysis date: 2026-07-12. This report supersedes the preliminary production summary in the
[implementation roadmap](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) and should be read with the
[paper-grounded expected results](../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/EXPECTED_RESULTS.md).

**Follow-up 2026-07-21:** the matched `S ≡ 0` Yang-geometry run completed cleanly with
`t½,0 = 60.2866`, passing the freshwater target inferred from Yang to about 2.2 %. Combined
with this production run, PARTIES gives `f̄5/f̄0 = 0.3352` versus Yang's plotted `≈0.5`.
The freshwater path is therefore validated directly; the unresolved Stage-A miss is now the
salt-specific normalized response or the published/digitized reference. See
[YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md).

**Allocation-constrained closure 2026-07-21:** the data-local matched-volume audit is now
complete for both `Sm=0` and `Sm=5`. It finds resolved salt recovery lengths (median
27.7--98.4 cells), no profiles thinner than four cells, leakage indices of only
`1.25e-4`--`1.92e-4`, and no heat-starved/refreezing rows. The slow salty melt follows reduced
heat flux and superheat, so the predeclared rule rejects the `Sc=500` branch. A subsequent
four-rank Le=100 equation gate tested Yang's finite-interface salt operator mapped onto PARTIES.
It was stable but worsened both the analytic liquid-salinity profile and Stefan-front constant
from both production-like and smooth nonzero-time starts. The refined gate and every production
branch are therefore stopped; see
[YANG_SALT_OPERATOR_GATE_REPORT.md](YANG_SALT_OPERATOR_GATE_REPORT.md).

**Route correction 2026-07-27.** The sentence that stood here — *"an accurate Yang reproduction
now requires the complete coupled Allen--Cahn/Hester phase--salt method, preferably in the
official AFiD-MuRPhFi code"* — is **withdrawn**. The goal is a PARTIES reproduction; an AFiD run
validates AFiD. No evidence identifies PARTIES' conservative Cahn--Hilliard formulation as the
defect, and against the one exact reference where the two salt formulations were compared head to
head, PARTIES won. The controlling plan is now
[YANG_ECCO_IMPLEMENTATION_ROADMAP.md](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) §A.5, which tests the
never-performed grid/interface-width convergence pair and the `f̄(Sm)` curve shape before any
interface-model rewrite.

## 1. Evidence and reference hierarchy

The benchmark definition comes from Yang et al., *Journal of Fluid Mechanics* 969 (2023), R2:

- [publisher article](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/ice-melting-in-salty-water-layering-and-nonmonotonic-dependence-on-the-mean-salinity/749F911A54CDD02671A7CA3597A4FB0A)
- [publisher PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/749F911A54CDD02671A7CA3597A4FB0A/S0022112023005827a.pdf/ice_melting_in_salty_water_layering_and_nonmonotonic_dependence_on_the_mean_salinity.pdf)
- [DOI 10.1017/jfm.2023.582](https://doi.org/10.1017/jfm.2023.582)
- [official AFiD-MuRPhFi repository](https://github.com/chowland/AFiD-MuRPhFi), for the reference
  phase-field/salinity solver

Local input and validation records are:

- [production input](../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/parties.inp),
  [job script](../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/jobscript.sh), and
  [expected results](../../PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_production/EXPECTED_RESULTS.md)
- [fresh Stefan gate](STAGE_A_STEFAN_GATE_REPORT.md),
  [salty Stefan gate](STAGE_A_SALTY_STEFAN_GATE_REPORT.md), and
  [EOS/ice-penalization gates](STAGE_A_EOS_ICE_PENALIZATION_GATE_REPORT.md)
- completed-run [time series](/bigscratch/mjalabert314/Yang_production/final_timeseries.csv),
  [summary figure](/bigscratch/mjalabert314/Yang_production/final_analysis.png), and
  [raw log](/bigscratch/mjalabert314/Yang_production/run.log)

The old `monitor_timeseries.csv` is not used for the verdict. Its volume calculation sliced HDF5
ghost cells incorrectly. The corrected `final_timeseries.csv` uses the physical domain and is
consistent with direct field integration.

## 2. Run completion and provenance

The production directory is `/bigscratch/mjalabert314/Yang_production`. The run used 144 MPI
ranks on the `1440 x 1440 x 1` `TWOD_CARTESIAN` grid, completed 400,019 iterations through
`t/t_ff = 200`, and took 246,733 s (2 d 20 h 32 min). It wrote 101 `Data_*` and 101 `Resume_*`
states spanning `t=0` through `t=200`. `production_job.err` is empty, and the raw log contains no
failed CH, velocity, concentration, or HYPRE solve.

The archived files are byte-identical to the repository copies used for this audit:

| Artifact | SHA-256 |
|---|---|
| `parties.inp` | `5345ad5f03ffd2fcea3a232af8b691cff0be5408aa16e923268623316e5eaf1e` |
| `jobscript.sh` | `cf664454b0dec47aea534ccf7992bf9a992671a86b97b2c0483556ed94389507` |
| `parties` executable | `36af49c0b517b22f1097aa14f0cabe5f8fdc1ebbeb0ab5ab0621df6d8727309a` |

This rules out an input-copy, executable-copy, early-stop, or missing-output explanation for the
benchmark discrepancy.

## 3. Input audit

The physical inputs match the intended Yang `Sm=5 g kg^-1`, `Delta Sv=5 g kg^-1` reference case.

| Item | Yang reference | Production value | Finding |
|---|---:|---:|---|
| Geometry | `Lx/H=1`; 3-D depth ratio `Ly/H=0.5` | 2-D unit square plus one storage slab | Correct for comparison with the paper's 2-D data |
| Ice | vertical block at `x>=0.9H`, thickness `0.1H` | `vof_slab_x0=0.9` | Correct |
| Initial temperature | water `20 degC`, ice `0 degC` | `theta_water=1`, `theta_ice=0` | Correct |
| Initial salinity | mean 5; top 2.5, bottom 7.5 `g kg^-1` | nondimensional top `0.5`, bottom `1.5`; zero in ice | Correct stable stratification |
| Walls | no slip, no heat flux, no salt flux | matching velocity and scalar boundary conditions | Correct |
| Thermal control | `RaT=10^7`, `Pr=10` | `Re=sqrt(RaT/Pr)=1000`, `Pe_T=10^4` | Correct |
| Solutal control | `Sc=1000`, `Le=100` | `Pe_S=10^6` | Correct |
| Stefan number | Yang `L/(cp Delta T)=4` | PARTIES reciprocal convention `stefan=0.25` | Correct |
| EOS | `Cb=0.011`, `b0=0.77`, `T0=4`, `cS=-0.25` | `betaT=1`, `betaS=1.75`, `Tmd0=0.2`, `Tmd_slope=-0.0625` | Correct nondimensionalization |
| Liquidus | `m=0.056 degC/(g kg^-1)` | `liquidus_slope=m Sm/Delta T=0.014` | Correct |
| End time | published volume-curve window through 200 | `time_max=200` | Correct |

The published text defines `Delta Sv = S_top-S_bot` while also stating `S_bot>=S_top` and plotting
positive `Delta Sv`. Those statements have an internal sign inconsistency. The production input
uses the physically and graphically consistent positive magnitude `S_bot-S_top=5`, with denser,
saltier water at the bottom. This is not a production-input error.

### Reference-grid correction

Yang's Figure 1 grid study is explicitly for `Sm=5 g kg^-1`, `Delta Sv=5 g kg^-1`: it selects
base `nx=288`, while the salinity/phase refined resolution is fixed at five times the base, hence
`1440`. The paper's `432 x 432` velocity/temperature and `2880 x 2880` salinity/phase grids apply
**specifically to `RaT=10^8`, `RaS=2x10^10`**, not this `RaT=10^7` reference case. The separate
3-D example uses `288 x 288 x 144` and `864 x 864 x 432`.

Therefore PARTIES' `1440^2` single grid matches the reference case's refined salinity/phase
resolution and is five times finer than its base velocity/temperature resolution. There is no
paper-supported 1440-versus-2880 resolution deficit for this run.

## 4. Acceptance results

| Required/health criterion | Acceptance | Result | Status |
|---|---:|---:|---|
| Half-melt time `t_1/2` | 106.2--129.8 (target about 118) | **179.8484** | **FAIL** |
| `V(200)/V0` | 0.27--0.37 (target about 0.32) | **0.4561566** | **FAIL** |
| Melt-front spacing | about `0.287H`, 3--4 layers, allowing layer ambiguity | `0.249--0.296H` in resolved interior layers; late upper merge `0.416H` | Qualified **PASS** |
| Top morphology | thicker/slower-melting ice at top | present | **PASS** |
| Salt conservation | relative drift at most `1e-8` | max `1.33e-15`; final `-1.11e-16` | **PASS** |
| Enthalpy | relative drift at most `1e-6` over diagnostic windows | max 2-unit-window drift `4.73e-8`; cumulative `-2.982e-6` over 200 | Window **PASS**; cumulative trend noted |
| Liquid fraction | exactly in `[0,1]` | exactly in `[0,1]` | **PASS** |
| Projection | divergence below `1e-6` | maximum `7.59e-10` over 1,200,054 projections | **PASS** |
| Solvers | converge without abort | all CH/velocity/scalar/HYPRE solves converged | **PASS** |
| Deep-ice velocity | strongly damped | at `t=200`: RMS `1.31e-4`, max `2.79e-4` | No instability; relevant to sensitivity study |

Because both physical comparisons are required, a morphology pass cannot compensate for the
failed global melt curve. The overall Stage-A decision is **FAIL**.

## 5. Ice-volume comparison and corrected slope interpretation

| `t/t_ff` | Digitized Yang `V/V0` | PARTIES `V/V0` | PARTIES minus Yang |
|---:|---:|---:|---:|
| 40 | 0.78218 | 0.82134 | +0.03916 |
| 60 | 0.69664 | 0.77427 | +0.07763 |
| 80 | 0.61855 | 0.72804 | +0.10949 |
| 100 | 0.55273 | 0.68015 | +0.12742 |
| 118 | 0.50360 | 0.63788 | +0.13428 |
| 120 | 0.49000 | 0.63327 | +0.14327 |
| 140 | 0.44218 | 0.58807 | +0.14589 |
| 160 | 0.39164 | 0.54298 | +0.15134 |
| 180 | 0.34818 | 0.49967 | +0.15149 |
| 198 | 0.31455 | 0.46049 | +0.14594 |

The discrepancy is **not a constant approximately 30 % slope deficit**:

- over `t=40--100`, Yang's digitized slope is `-0.003827`, whereas PARTIES gives
  `-0.002331`, about 39 % too slow;
- over `t=120--200`, Yang gives `-0.002292` and PARTIES gives `-0.002210`, which are nearly equal.

Most of the permanent volume offset is created during `t=40--120`, then carried into the late
window. This is why a sensitivity restart at `t=100` would begin too late to diagnose the cause.
The suitable common restart is `Resume_10.h5` at `t=20`.

Using half-melt rates, `f_PARTIES/f_Yang = 118/179.8484 = 0.6561`: the production case is 34.39 %
slow in the benchmark's mean-rate metric.

## 6. Morphology result

The `F=0.5` contour gives representative adjacent-scallop spacings:

| Time | Measured spacings (`H`) | Interpretation |
|---:|---:|---|
| 52 | 0.296, 0.278 | close to the `0.287H` prediction |
| 100 | 0.284, 0.290 | quantitative agreement |
| 146 | 0.270, 0.263 | modest thinning as the field evolves |
| 200 | 0.249--0.256 interior; 0.416 upper merged interval | interior layers persist; upper layers have merged |

The simulation also develops the expected thicker/slower-melting upper ice associated with the
fresh layer under the lid. This establishes a qualified morphology pass: the double-diffusive
layering mechanism is present and quantitatively close before late layer merging.

## 7. Why the benchmark did not pass

The evidence excludes several simple explanations. The input parameters and reference grid are
correct; the run completed; budgets and solvers are healthy; the final physical-domain integration
is correct; and reconstructed `dF/dt` agrees with the stored melt source to about 0.2--1 %. The
melt-source normalization is approximately 4.8 % stronger locally and 11.1 % stronger after
integration than the corresponding idealized reference normalization, so an accidentally weak
source is not supported and would in any case have the wrong sign to explain the slow melt.

The remaining causes should be tested in this order:

1. **Leading hypothesis: ice penalization is outside the validated production regime and differs
   materially from Yang.** Yang uses `eta=dt` in `-phi u/eta` and directly forces velocity to zero
   for solid fraction above 0.9. PARTIES uses a fixed implicit `darcy_tau=1e-4` throughout the
   diffuse ice fraction and has no cutoff. With the production step near `5e-4`, the nominal
   damping coefficient is about five times stronger wherever the same phase weight is compared.
   More importantly, the production Brinkman thickness is only about `delta/dx=0.455`, whereas the
   penalization gate covered `delta/dx={2,4,8}`. Thus this case extrapolates beyond that gate and can
   suppress the melt-feeding interfacial circulation during precisely the `t=40--120` window.
2. **Secondary hypothesis: phase-field/melt-coupling difference.** Yang evolves an Allen--Cahn
   phase field with its Hester coupling; PARTIES uses its validated conservative Cahn--Hilliard
   transport and superheat/source mapping. The preliminary Stefan gates validate isolated fronts,
   but do not prove identical convective interface response in the full double-diffusive problem.
3. **Not the leading hypothesis: spatial resolution.** The corrected paper reading shows that
   `1440^2` is the exact case's refined resolution, not half of a `2880^2` target. A resolution
   study may still quantify PARTIES discretization uncertainty, but it should not precede the
   model-form sensitivity tests.

The cumulative enthalpy trend of `-2.982e-6` is above the full-run `1e-6` aspiration but is two
orders of magnitude too small to explain a roughly 0.14 absolute ice-volume offset at `t=120`.

## 8. Solutions and next steps

> **Current controlling plan (2026-07-21, supersedes the exploratory sequence below): STOP all
> parameter and salt-only branches.** The paired transport audit rejected `Sc=500`, and the
> Le=100 salt-operator gate rejected both refinement and a PARTIES production test: the mapped
> Yang salt equation doubled the liquid-profile error from a smooth analytic start and increased
> the front-constant error. Do not repeat Darcy, CH-mobility, `Sc`, interface-thickness, or source
> prefactor sensitivities. For an accurate Yang reproduction, first reproduce the official
> AFiD-MuRPhFi 1-D multicomponent validation at coarse and refined resolution, then run only one
> checkpoint-extendable 2-D case to `t=60`; continue that same trajectory to `t=200` only if it
> closes at least 25 % of the present `t=60` volume gap with healthy conservation. If Yang must be
> reproduced inside PARTIES, port and validate the complete Allen--Cahn/Hester phase--salt update
> before any further physics run. Exact results and gates:
> [YANG_SALT_OPERATOR_GATE_REPORT.md](YANG_SALT_OPERATOR_GATE_REPORT.md).

### Step 1 -- three controlled penalization restarts

Branch from the common `t=20` state, `/bigscratch/mjalabert314/Yang_production/Resume_10.h5`, and
run to `t=80` with all fields, executable, and inputs unchanged except restart/end-time settings
and `darcy_tau`:

| Branch | `darcy_tau` | Purpose |
|---|---:|---|
| control | `1.0e-4` | proves restart continuity against the original trajectory |
| intermediate | `2.5e-4` | measures monotonic sensitivity |
| Yang-timescale proxy | `5.0e-4` | matches the dominant production `dt` scale without yet adding Yang's cutoff |

Set `resume=1`, link/copy the common checkpoint as the branch's restart input, retain the same
144 ranks and `1440^2` grid, set `time_max=80`, and save at least every two free-fall units.
Compare `V/V0` and interval slopes over `t=40--80`, contour spacing, near-front velocity, deep-ice
velocity, salt, enthalpy, and source-reconstructed `dF/dt`. Do not use the old ghost-sliced monitor.

#### Submitted branch record (2026-07-12)

The isolated branches are under
`/bigscratch/mjalabert314/Yang_penalty_branches_t20/`. Each contains a private copy of
`Data_10.h5`, a private copy of `Resume_10.h5` renamed to `Resume.h5`, and the exact production
executable. The completed `/bigscratch/mjalabert314/Yang_production` directory was not modified;
its `Resume.h5` still points to `Resume_100.h5`.

| Branch directory | `darcy_tau` | SLURM job | State at submission |
|---|---:|---:|---|
| `tau_1e-4` | `1.0e-4` | `277955` | pending (priority) |
| `tau_2p5e-4` | `2.5e-4` | `277956` | pending (priority) |
| `tau_5e-4` | `5.0e-4` | `277957` | pending (priority) |
| dependent analysis | -- | `277958` | pending on jobs 277955--277957 |

All three inputs have `resume=1` and `time_max=80`. The branch seed records `time=20`,
`dt=5e-4`, and `noutput=10`, so the first new files will be `Data_11.h5` and `Resume_11.h5`.
Checksums of the private `Data_10.h5`, `Resume.h5`, and `parties` copies match the originals.
When all three simulation jobs complete successfully, job `277958` will run the preserved
`analyze_yang_production.py` workflow and write the comparison CSV files under the branch root's
`analysis/` directory.

#### Branch outcome (2026-07-13): darcy_tau is NOT the cause

Interim analysis of the two completed branches plus the partial third
(`analysis_partial/` under the branch root; `tau_branches_overlay.png`), run with
`analyze_yang_production.py` against the production baseline:

| Quantity | control `1e-4` | `2.5e-4` | `5e-4` (partial, t→24) | Yang digitized |
|---|---:|---:|---:|---:|
| `V/V0` at t=40 | 0.82134 | 0.82153 | -- | 0.78218 |
| `V/V0` at t=60 | 0.77427 | 0.77421 | -- | 0.69664 |
| `V/V0` at t=80 | 0.72804 | 0.72717 | -- | 0.61855 |
| slope t=40--80 | -0.002332 | -0.002359 | -- | ~-0.0038 |
| max solid speed | 5.07e-4 | 1.23e-3 | 2.34e-3 | -- |
| salt drift | 9.9e-16 | 9.9e-16 | 3.7e-16 | -- |

- The control branch reproduces the production trajectory to `1.4e-8` in `V/V0` at `t=80`
  (identical to all printed digits at t=40/60): restart continuity is proven.
- `darcy_tau=2.5e-4` closes **0.8 %** of the t=80 volume gap (0.00087 of 0.10949); its
  t=40--80 slope steepens by 1.1 %. This is within restart noise.
- The partial `darcy_tau=5e-4` branch (Yang `eta=dt` timescale proxy) overlays the control
  through t=24 (`+9e-5` in `V/V0`, the wrong sign to close the gap).
- Solid creep grows in proportion to `tau` (5.1e-4 -> 2.3e-3) while melt does not change:
  the penalization is genuinely weakened five-fold, and the melt-feeding circulation does
  not respond. The `delta/dx=0.455` extrapolation concern is therefore answered empirically:
  the Brinkman drag is not throttling the interfacial circulation.

**Decision per Step 2, third bullet:** keep `darcy_tau=1e-4`; the leading hypothesis moves
to the Allen--Cahn (Hester coupling, `D=1.2*kappa`, `eps=dx`, `eta=dt`, forcing cutoff at 0.9)
versus Cahn--Hilliard (`Cn=0.75dx`, `Pe_CH=0.9/Cn`, superheat source `(St/Pe_T)(theta-theta_L)/eps`
with `w(F)=F(1-F)/(sqrt(2)Cn)` deposition) interface-coupling difference under convection,
which the quiescent 1-D gates by construction could not discriminate.

#### CH-mobility branch outcome (2026-07-16): mobility is NOT the cause either

Two `Pe_CH` branches from the same t=20 seed (`/bigscratch/mjalabert314/Yang_mobility_branches_t20`,
record + CSVs in `analysis/`): `pech_8640` (5x slower band relaxation, matching Yang's
`D = 1.2*kappa_T` diffusivity scale within 4%) ran to t=52 and tracks the control to
`+1.4e-4` in `V/V0` at t=40 -- through the steepest part of the deficit window -- with a
fully healthy band (F exactly in [0,1], CH iterations <= 6, salt drift 9e-16). The 2x
intermediate (`pech_3456`) was made redundant and both jobs were cancelled at that point to
free the queue.

**Interface machinery exonerated (2026-07-16, decisive).** Direct measurement on the
production snapshots shows the CH interface is quantitatively exact at production parameters:
(i) the |grad F|-weighted band superheat equals the kernel requirement `V_n*Pe_T*eps/St` at
every sampled time (0.0125 vs 0.0143 at t=4; 0.0044 vs 0.0050 at t=100) -- the interface sits
at the liquidus to ~0.5-1% of Delta-T; (ii) the fitted thermal diffusivity from the front
profile is `alpha_eff/alpha_in = 1.15-1.16` at t=2-4 (excess = incipient convection), profile
shape erfc-exact (`early_profile_audit.png`); (iii) the early quasi-diffusive front law matches
the correct **reservoir-driven** similarity solution `sqrt(pi)*lam*exp(lam^2)*(1+erf(lam)) = St`
=> `lam_exact = 0.122` for `stefan = 0.25`: measured `lam_inst = 0.127` at t=2 (+4%, rising
with convection). (Note: the hot-wall gate solution `lam = 0.3401` does NOT apply to the
reservoir-driven production geometry -- an earlier comparison against it was retracted.)
Combined with the tau and Pe_CH nulls: **the deficit is not in the phase-change/interface
machinery.** Remaining causes: convective heat delivery to the front (double-diffusive
film/roll structure), or bias in the digitized reference itself. Discriminator chosen: the
salt-free melting Rayleigh-Benard benchmark (Favier, Purseed & Duchemin, JFM 858, 2019) run
with the production executable -- see `ECCO_TESTS/StageA_Yang/MeltingRB/`.

### Step 2 -- choose the model correction from the sensitivity result

- If increasing `darcy_tau` moves the curve monotonically toward Yang without unacceptable ice
  motion, extend the best branch to `t=200` and reapply both primary acceptance tests.
- If the larger `darcy_tau` improves melt rate but loses rigidity, implement and validate Yang's
  solid-fraction direct-forcing cutoff before repeating the branch. Re-run a penalization gate at
  `delta/dx<1`, because the existing `{2,4,8}` sweep does not cover production.
- If the three curves are insensitive to `darcy_tau`, keep penalization fixed and isolate the
  Allen--Cahn versus Cahn--Hilliard/interface-coupling difference with a short matched case. Audit
  the phase-field thickness, mobility, liquidus pinning, and convective source mapping together;
  changing a source prefactor alone is not justified by the present source-reconstruction result.

### Step 3 -- convergence and closure

After selecting the model form, perform one spatial/interface-thickness convergence comparison.
The relevant baseline is Yang's `288/1440` multiple-resolution case, not `432/2880`. Then rerun
through `t=200`. Stage A closes only when both `t_1/2`/`V(200)/V0` and the layer-spacing criterion
pass while conservation and solver-health checks remain satisfied. Until then, do not advance to
the Stage-B particle/ECCO production campaign.
