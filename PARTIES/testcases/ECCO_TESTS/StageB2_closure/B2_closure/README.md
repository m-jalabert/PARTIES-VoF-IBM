# B.2 sediment-from-ice validation

The controlling result is [the closure report](../../../../../docs/md_files/NASA-JPL/YANG_ECCO_B2_CLOSURE_REPORT.md).
These are synthetic diagnostics. They do not specify or authorize the real ECCO scenario.

## Build and configuration

The tested changes are in the shared `PARTIES/src` tree. Both new switches,
`VOF_DIFFUSE_SEDIMENT_ICE_TRANSPORT` and `ECCO_PROFILES`, are **off by default**.
`Boundary.validation.h` records the exact sediment configuration used by the gates:
3-D uniform Cartesian grid, x/z periodic, y walls, resolved IBM, symmetric sediment
CH mask, phase change, three scalars, hard-threshold ice penalization, and profiles.
`Boundary.default_before.h` preserves the pre-promotion default configuration.
`solver_changes.patch` records just this task's solver delta against the coworker's
working tree, excluding the subsequently added default-off switch declarations.

On Anvil, from this directory:

```bash
module load gcc/11.2.0 openmpi/4.1.6 hdf5/1.10.7 fftw/3.3.8
python build_validation.py --configuration sediment --destination build_final_new
```

The build script rejects an existing destination, clones the shared source,
applies the recorded configuration, and writes a source/binary SHA-256 manifest.
It does not run a simulation. `final_build_manifest.json` records the clean
rebuild that was checked against the original candidate: all saved seed fields
were bit-identical despite differing executable hashes.
`--configuration remap` builds the manufactured translation fixture by injecting
`remap_fixture.inc` into the isolated copy only. Never use that test executable
for a physical simulation. `audit_campaign.py` refreshes `campaign_jobs.json`
from Slurm accounting, including rejected and cancelled trials.

All simulation output belongs under `/anvil/scratch/x-mjalabert/ECCO_StageB/`.
Each job script creates a unique `B2closure_<arm>_<jobid>` directory and captures
the input files, executable, and executable hash there. Submit **from the arm's
directory**, because the scripts read `SLURM_SUBMIT_DIR`:

```bash
cd coupled16
sbatch job.sh
```

The archived scripts name the immutable candidate executables used for their
recorded jobs. For a new rebuild, set that executable path to the newly built
binary before submission. Do not reuse an output directory or an incompatible
checkpoint. The final transport checkpoint version is **2**; legacy liquid-CH
and version-1 periodic-mask checkpoints are deliberately rejected.

## Accepted validation inputs

| Input directory | Purpose |
|---|---|
| `coupled16` | Complete lock/release/settle/collision at d/dx=16, tau=0.001, nonzero passive salt |
| `coupled16_halfdt` | Half maximum timestep, doubled ramp steps to preserve the initial release-ramp duration |
| `coupled_v9` / `coupled_v9_mpi8` | Cheap d/dx=8 coupled MPI comparison, including periodic particle wrap |
| `settle16`, `settle16_nodarcy`, `settle16_reference` | Genuine ice-free settling: stiff Darcy, negligible Darcy, conventional IBM controls |
| `settle24`, `settle32` | Transient settling refinement through t=6 |
| `hot_water` | Warm water with phase change enabled must not melt the rock or create ice/tracer |
| `seed_v9`, `restart_v9_*` | Closed-box nonzero-salt budgets and exact restart reproduction |
| `remap8`, `remap16` | Manufactured periodic remap; uses the explicitly test-only `remap_fixture.inc` injection |
| `stagea_final` | Final default-off Stage-A regression |
| `restart_guards` | Expected MPI aborts for wrong transport version and missing Data file |
| `coupled16_halfdt`, `coupled16_quarterdt` | Matched-ramp timestep triple, **old** ramp |
| `coupled16_rampfix{,_halfdt,_quarterdt}` | Same triple under the release-ramp fix |
| `coupled16_mpi32` | 32-rank twin of `coupled16`: roundoff-perturbation floor |
| `coupled_v9_rampfix{,_mpi8}` | MPI-agreement gate re-run under the fix |
| `settle24_dt10`, `settle24_dt05` | **Ice-free settling timestep control** (positive control) |
| `coupled16_smooth{,_halfdt}` | Hard ice threshold DISABLED. Diagnostic only, see below |

`coupled16_smooth*` uses `--configuration smooth`, which injects
`Boundary.smooth_penalization.h` (the validated header with
`VOF_DIFFUSE_ICE_PENAL_THRESHOLD` commented out). It exists solely to measure what
that discontinuity costs. **It is not a validated configuration**: it delays
release by 10.3 time units and worsens the ice-rigidity gate. Never use it for a
physical scenario.

Older `tiny_*`, `coupled`, `coupled_v6*`, and `remap_replay` directories retain
rejected experiments. Their output is evidence of failed hypotheses, not an
alternative accepted implementation. `coupled_v9_halfdt` changed the physical
duration of the step-counted release ramp and is a sensitivity diagnostic, not
a clean temporal convergence comparison.

## Reproduce the audits

Python dependencies: NumPy, h5py; Matplotlib only for figures. No MPI or SU charge
is needed to analyze saved fields.

```bash
python analyze_profiles.py /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled16_20478206 --output coupled16_profiles
python verify_outputs.py profiles /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_coupled_v9_20477868
python verify_outputs.py remap /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_remap8_20478180 /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_remap16_20478181
python verify_outputs.py restart /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_restart_v9_continuous_20478188 /anvil/scratch/x-mjalabert/ECCO_StageB/B2closure_restart_v9_resume_20478192

# B.2 gate table (checks 1, 2, 4).  --check validates the reimplementation
# against the recorded summary before it is trusted on a new run.
python analyze_motion.py <run> --name coupled16 --check motion_summary.json

# Timestep refinement.  Separates the ill-conditioned escape lag from the
# well-conditioned open-water settling; three runs give an observed order.
python refine_timesteps.py <dt0.01> <dt0.005> <dt0.0025> --dt 0.01 0.005 0.0025

# Confirm a solver change altered nothing before the grain was released.
python verify_prerelease.py <old-run> <new-run>
```

`analyze_motion.py` re-derives the gate table that was previously generated inline
and not preserved; it reproduces `motion_summary.json` to 1.1e-16.

The audits exclude duplicated high boundary planes. Profiles distinguish
liquid-weighted means/fluxes from whole-mesh conserved scalar integrals.
The independent dissipation check covers interior y rows and all periodic x/z
cells; it does not independently verify the wall derivative stencil.
Budget integration currently supports the tested unit thermal diffusivity ratio.
The scalar values within the IBM's excluded sediment support are extension
values; physical tracer extrema must be evaluated outside that support.

**Timestep validity.** Ice-free settling is converged to 0.28% in `U_y` at
`max_dt = 0.01` (`settle24_dt10/dt05`). The grain's *escape from the ice* is not
converged at any affordable timestep and its error is a one-signed bias, so
release time, escape duration and release-to-contact timing are not usable
quantities. See the closure report before using any release-relative number.

The short settling column establishes transient refinement and wall collision,
not a terminal-velocity plateau or an unbounded drag correction. Scenario-specific
EOS, mesh, timestep, and interface-width choices require the selected physical
scenario. No production mixing efficiency is inferred from the synthetic BPE.
