# Stage A nonlinear-EOS and ICE-penalization gate report

**Status:** PASSED after periodic buoyancy endpoint correction, 2026-07-05.

## Scope

This report closes the two steady-flow Stage-A gates that complement the
quiescent pure and salty Stefan tests:

- `EOS_NONLINEAR`: reconstruct the imposed Roquet buoyancy pointwise from the
  curvature of a steady unidirectional velocity profile.
- `ICE_PENALIZATION`: compare steady Brinkman-Poiseuille profiles against the
  discrete diffuse-interface BVP for Darcy layers of 2, 4, and 8 grid cells.

The working cases are `PARTIES/testcases/1DNonlinearEOS/` and
`PARTIES/testcases/1DIcePenalization/`.

## Defect isolated before certification

The first completed runs consistently produced 0.75 of the nominal velocity
amplitude. An independent all-liquid fixed-pressure-gradient Poiseuille test
returned unit amplitude, excluding the momentum/viscous operator as the cause.

Inspection of `Velocity_add_buoyancy_2_RHS` then showed that its component loops
did not include the periodic wrap face. Both gates use gravity along periodic
`z` with `NZM=4`, so only three of four physical `w` faces received buoyancy.
Projection produced the observed factor `3/4`.

The correction gives the `u`, `v`, and `w` buoyancy loops component-specific
periodic endpoints matching the normal momentum RHS. Cell-centred EOS assembly
and non-periodic loop bounds are unchanged.

## EOS gate

The original converged states at `t=12` were resumed with the corrected
executable and relaxed to `t=20`.

| Metric | Thermal branch | Salty branch |
|---|---:|---:|
| Original measured/BVP amplitude | 0.74999418 | 0.74999768 |
| Corrected measured/BVP amplitude | **0.99986457** | **0.99991996** |
| Nominal BVP relative error | 1.41e-4 | 8.39e-5 |
| Scaled profile residual | 2.59e-5 | 6.26e-5 |
| Last-output relative change | 6.05e-5 | 3.64e-5 |
| Median reconstructed/input buoyancy | 0.99991565 | 1.00000228 |

The reconstructed thermal parameters are:

- `betaT = 1.0002443` (nominal 1),
- `Tmd0 = 0.4000146` (nominal 0.4),
- `q = 2.0004528` (nominal 2).

The reconstructed salty polynomial coefficients are:

- `c0 = -0.3600050` (nominal -0.36),
- `c1 = 1.1499418` (nominal 1.15),
- `c2 = -0.2499509` (nominal -0.25).

**EOS conclusion:** the nonlinear thermal anomaly, density-maximum location,
quadratic exponent, salinity term, and salinity shift all pass both functional
and absolute-amplitude checks.

## ICE-penalization full sweep

The original `t=4` states were resumed with the corrected executable and
relaxed to `t=6.4` using `dt=0.002`. The BVP uses the measured diffuse `F` and
salinity fields and the correct gravity sign:

`nu*w'' - (1-F)/darcy_tau*w = 2*s`.

| Metric | 2 dx | 4 dx | 8 dx |
|---|---:|---:|---:|
| Original measured/BVP amplitude | 0.75002078 | 0.75000453 | 0.75000086 |
| Corrected measured/BVP amplitude | **0.99999698** | **0.99997370** | **0.99994732** |
| Nominal BVP relative error | 3.74e-6 | 2.72e-5 | 5.45e-5 |
| Scaled profile residual | 1.28e-6 | 1.36e-6 | 2.62e-6 |
| Last-output relative change | 1.33e-5 | 3.22e-5 | 5.65e-5 |
| Discrete decay length | 0.00789248 | 0.01566551 | 0.03127032 |
| Fitted decay length | 0.00789271 | 0.01566558 | 0.03127038 |
| Decay error vs discrete theory | 2.91e-5 | 4.06e-6 | 1.89e-6 |
| Far-ice leakage | 4.19e-7 | 2.52e-5 | 2.59e-4 |

Across the sweep, the measured interface-slip exponent is `0.59727023`; the
diffuse-BVP exponent is `0.59728621`, an absolute difference of `1.60e-5`.
The exponent differs from the ideal sharp-interface value because the chosen
diffuse band is not asymptotically thin relative to the 2-dx Darcy layer; the
measured data follow the appropriate diffuse discrete reference.

A corrected-source rank-count check resumed the same 4-dx state from `t=4` to
`t=4.2` on one and four ranks. The worst stored-field difference was
`1.92e-12` (pressure); the `w` difference was `1.29e-13`.

**ICE-penalization conclusion:** all three cases pass absolute amplitude,
profile, steady-state, decay-length, leakage, and sweep-scaling checks. The
implicit Darcy coefficient and its viscosity-relative normalization are
correct.

## Reproducible artifacts

Combined periodic-fix validation:

- `PARTIES/testcases/periodic_fix_validation/analyze_periodic_fix.py`
- `periodic_fix_results.csv`, `periodic_fix_profiles.csv`, and comparison plots
  in that directory.

Corrected complete penalization sweep:

- `PARTIES/testcases/1DIcePenalization/analyze_corrected_sweep.py`
- `corrected_penalization_results.csv`
- `corrected_penalization_profiles.csv`
- `corrected_penalization_history.csv`
- `corrected_penalization_profiles.png`
- `corrected_penalization_slip_scaling.png`
- `corrected_penalization_history.png`

## Stage-A decision

The four quantitative module gates now pass: pure Stefan, salty Stefan,
nonlinear EOS, and ICE penalization. Before committing resources to the full
1440 x 1440 x 4 Yang production run, the remaining roadmap prerequisite is the
small full-flag integrated A.2 shakeout. It should combine phase change,
temperature, salinity, nonlinear buoyancy, and ice penalization in the actual
Yang orientation and verify stability, conservation, salt exclusion, and MPI
behaviour. If that integrated shakeout passes, proceed directly to A.4 Yang.
