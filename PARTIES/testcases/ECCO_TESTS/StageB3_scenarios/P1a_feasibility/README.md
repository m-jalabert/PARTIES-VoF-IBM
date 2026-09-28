# P1 idealized scenario — bounded feasibility segment, 2026-09-25

> **Superseded 2026-09-28 by [P1b](../P1b_shelf/README.md)**: the net freezing here is a diffuse-band
> artifact at liquidus slope ≈ 1, and melt-out is unaffordable; see
> [YANG_ECCO_B3_P1b_SHELF.md](../../../../../docs/md_files/NASA-JPL/YANG_ECCO_B3_P1b_SHELF.md) §1.


**Results audited 2026-09-27:** job 20906327 completed for **182.76 SU**.
Conservation passes; no release occurred, the model has slight net freezing,
and the relative rigidity criterion fails. Full P1 remains open. See
[the result report](RESULTS_20906327.md), [diagnostics](results_20906327/P1_diagnostics.png),
and the compact machine-readable archive in `results_20906327/`.
The original submission plan below is retained as provenance.

The user authorized P1 and accepted the physical case proposed on September 25:
1 mm mineral grain, 2650 kg/m³, water/bottom CT=0 °C, interface SA=34 g/kg,
stable gradient 0.01 (g/kg)/m, p=0 dbar. This is an idealized near-surface
millimetre column, not a calibrated ocean site. The very small ambient salinity
change across P1 (0.00008 g/kg below the interface) is intentional; fresh melt
provides much stronger local stratification. Ice starts at ambient freezing.

`prepare.py` reproduces the numeric inputs and `case.json`. All density fits use
GSW exact Gibbs density with explicit CT-to-in-situ conversion. The broad fit
includes dilute and diffuse-band states, even those below local freezing; it
does not assert thermodynamic stability of supercooled pure liquid there.
Physical material constants are declared idealizations. CT uses cp0; latent
heat is fixed. Sc=70 gives a modeled salt diffusivity 25.7 times the declared
1e-9 m²/s comparison value. Ice/water diffusivity ratio remains one.

## EOS decision before simulation

The roadmap's linear starting candidate is rejected: a global linear density
fit over fresh-to-salty water gives positive d(rho)/dCT, incorrect in ambient
seawater. A density-only nonlinear least-squares fit can also conceal large
thermal-response errors. Instead, the existing q=2 EOS is selected, with its
density-maximum locus fitted to GSW, curvature fitted to d(rho)/dCT, and then
haline amplitude and additive constant fitted to density. There is no solver
change. The constant is absorbed by pressure. The haline velocity scale uses
the explicit fitted haline coefficient; the quadratic contributes a small
additional local salinity derivative. The usual linear Richardson mapping
is therefore inactive; `richardson={0,0,0}` is intentional.

`case.json` records holdout errors, local derivative errors, the freezing secant
residual, and all nondimensional coefficients. This is a checked approximate
EOS, not exact TEOS-10 in the flow solver. It preserves opposite thermal
buoyancy directions for cold freshwater and ambient seawater. Local thermal
errors near the density maximum must be read on an absolute scale.
The liquidus secant passes exactly through SA=0 and SA=34; its nonzero interior
error is part of the approximation and must accompany physical interpretation.
Retain this same EOS for any later P2/production, unless the campaign is
explicitly redesigned before those runs.

## Resource and observation decision

Allocation before submission: **71,699.3 SU**, no active jobs. P1 uses the planned
96×240×96 mesh and 4×10×4 domain. The physical case has Re≈8.96, PeT≈115.26,
St≈0.02213, G*=37.68. B.2's fast synthetic melt-out does not transfer.
The initial sensible heat can melt only 0.177d of a planar layer. A heat-only
estimate for 1d retreat, spending that reservoir first and using an 8.5d water
conduction path, gives ~36,435 nondimensional time units (~38 minutes physical).
At dt=0.002, the historical throughput anchors imply 0.36–2.42 million SU,
before contingency. This is an order-of-magnitude estimate, not a release-time
prediction or rigorous bound: cold-top heat loss, salinity diffusion and
convection change it. It is sufficient reason not to commit a full melt-out run.

The authorized P1 starts as a **bounded feasibility segment**: t≤10 (~0.62 s),
128 cores, at most 3 hours/**384 SU**, normal stop/checkpoint requested at 170
minutes. The nominal 5000 steps cost 98–663 SU; with 30% contingency, 128–863 SU.
The hard cap deliberately stops short if the slow anchor applies. Completion
of this segment is not completion of P1's post-escape objective. It tests the
actual salt/EOS/liquidus startup, initial geometry, budgets, melt direction and
throughput. There is **no automatic continuation**, P2, or production submission.
No quantitative release clock or terminal-velocity claim is permitted.

dt_start=0.0002, dt_cap=0.002; the source scale is ~0.9005. A planning speed of
5 gives advective CFL=0.24, below 0.3. The adaptive CFL/contact restrictions
remain active. Fifty ramp steps target 0.1 time units only if the cap is reached;
actual ramp duration must be measured if release occurs. Field cadence 0.025
limits displacement to 0.125d at that planning speed. About 250 GiB is reserved
as a conservative full-segment output estimate. All output and the isolated
source/build live in scratch; no shared source or executable is modified.

For a matched t=10 P2, 10,000 nominal steps at dt=0.001 cost 465–3145 SU
(605–4089 with contingency); this short window would still not establish the
post-escape comparison. No P2 numeric deck is represented as accepted.

## Provenance and checks

Build with the existing B2 `build_validation.py --configuration sediment`, to
`/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_build_20260925`. The exact header is
`Boundary.scenario.h` (identical to the validated sediment header, including
EOS_NONLINEAR). The promoted release-increment ramp remains present. Build
manifest, binary hash, complete flat inputs, case record, audit scripts and
Slurm metadata are copied into the unique job scratch directory.

Run `python prepare.py` to reproduce inputs; submit `sbatch job.sh` from a staged copy of this
directory under `ECCO_StageB3`. Do not re-submit to continue: inspect results, allocation and all
checkpoint/history files first. Heat auditing uses B2's full-physical-mesh
normalization and unit-diffusivity boundary-flux integration. Physical scalar
extrema exclude sediment support (`C_S<0.05`).

References: [GSW density](https://teos-10.github.io/GSW-Python/density.html),
[GSW freezing and temperature functions](https://teos-10.github.io/GSW-Python/gsw_flat.html).

Submission: **20906327**, initially PENDING (Priority). No runtime acceptance
is claimed. To audit it after complete output blocks exist:

```bash
python audit.py /anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_20260925_20906327 --output results.json
```

The audit reader was smoke-checked against archived B2 job 20492853; that is
reader verification only. The result is in the isolated build directory.

Storage policy (user instruction, 2026-09-25): all future B.3 builds, submission
files, scheduler logs and simulation output belong under
`/anvil/scratch/x-mjalabert/ECCO_StageB3/`. Historical B.2 runs stay in
`ECCO_StageB`. The build directory was moved without changing its binary.
The live Slurm configuration has no `short` partition; `part-debug` limits
wall time to 2 hours (and 2 nodes / 256 CPUs, one job per user), so this
3-hour segment retains `shared`. Physics, timestep, duration and SU cap are unchanged.

Replacement history: original job **20906230** was held and cancelled while
still pending; replacement **20906327** was submitted from
`/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/B3_P1_submit_20260925`.
