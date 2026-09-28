# Reduced integrated Yang checkout

**Status:** PASSED on 2026-07-05.

## Scope

This is the Stage-A A.2 code shakeout, not the full Yang physics benchmark.  It
uses the complete Yang compile-time flag set, two scalar fields, nonlinear EOS,
phase change, phase-weighted diffusion, and implicit ice penalization in the
Yang orientation (`x/y` no-slip, thin periodic `z`).

The final SLURM matrix was job `273214` and completed with exit code `0:0` in
102 seconds on one node (8 allocated tasks).  Raw output is in
`/bigscratch/mjalabert314/Yang_checkout_complete`.

| Run | Grid | MPI ranks | Purpose |
|---|---:|---:|---|
| `run_nz4_1rank` | 128 x 128 x 4 | 1 | coupled reference |
| `run_nz4_4rank` | 128 x 128 x 4 | 4 | MPI consistency |
| `run_nz8_4rank` | 128 x 128 x 8 | 4 | periodic-depth-count consistency |
| `run_static_control` | 128 x 128 x 4 | 4 | melting off, gravity off, negligible Darcy damping |

Each run advances 50 fixed steps at `dt=1e-3` to `t=0.05`.

## Results

All runs completed without NaN/Inf, convergence errors, or out-of-bound liquid
fraction.  The maximum accepted velocity-divergence residual was `9.97e-7`
against a `1e-6` projection tolerance.  Velocity and CH solves required at most
7 and 8 iterations, respectively.

| Check | Result | Gate | Status |
|---|---:|---:|---|
| 1-rank vs 4-rank worst stored-field difference | 1.339e-13 | < 1e-9 | PASS |
| NZ=4 vs NZ=8 worst z-mean difference | 2.001e-9 | < 1e-7 | PASS |
| Coupled-run maximum salt relative drift | 3.701e-16 | < 1e-8 | PASS |
| Coupled-run maximum enthalpy relative drift | 1.458e-8 | < 1e-6 | PASS |
| Coupled-run maximum deep-ice salt change | 3.069e-7 | < 1e-6 | PASS |
| Coupled liquid-volume increase | 1.065e-4 | > 0 | PASS |
| Static-control salt drift | 3.701e-16 | < 1e-12 | PASS |
| Static-control heat drift | 1.234e-16 | < 1e-12 | PASS |
| Static-control maximum speed | 0 | < 1e-12 | PASS |
| Static-control front displacement | 9.187e-6 (0.00118 cells) | < 1e-4 | PASS |
| Static-control liquid-volume drift | 1.445e-8 | < 1e-7 | PASS |

The static control has a local `max|F(t)-F(0)|=0.00510` from equilibration of
the analytic tanh profile under Cahn-Hilliard dynamics.  Its front and global
volume remain effectively stationary; this local relaxation is not melting or
advection.

## Checkout finding and correction

The first launch exposed an off-by-one pressure-loop abort: the eleventh
projection reduced divergence below `1e-6`, but the code aborted immediately
afterward because it checked only the iteration count.  `Temporal_int.c` now
aborts after that correction only when the residual is still above tolerance.
The corrected matrix completes normally.

## Decision

The reduced integrated A.2 checkout passes.  It demonstrates that the already
validated EOS and ice-penalization implementations operate together with phase
change, salinity, MPI decomposition, and thin-periodic geometry without a new
integration failure.  It does not validate Yang melt-front or ice-volume
physics; those are the acceptance criteria for the full A.4 production run.

## Reproducible artifacts

- `parties_nz4.inp`, `parties_nz8.inp`, `parties_control.inp`
- `jobscript.sh`
- `analyze_yang_checkout.py`
- `yang_checkout_summary.csv`
- `yang_checkout_timeseries.csv`
- `yang_checkout_interface.csv`
- `yang_checkout_fields.png`
- `yang_checkout_timeseries.png`
- `yang_checkout_consistency.png`

