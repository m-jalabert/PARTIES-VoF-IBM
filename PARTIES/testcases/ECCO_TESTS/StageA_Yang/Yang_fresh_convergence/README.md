# Yang freshwater convergence ladder (roadmap §A.5.2)

**Question.** PARTIES melts the Yang freshwater case (`Sm = 0`) in `t½ = 60.287`; Yang's
directly digitized 2-D curve gives `t½ = 107.42` — **1.78× too fast**. Is that a resolution /
interface-width effect, or a model (or reference) discrepancy?

PARTIES has never run a convergence study of this case. The completed `Sm=0` run inherited
`1440²` from the salty production input, where the mesh was set by `Sc = 1000`. With `S ≡ 0`
that constraint does not exist — only `Pe_T = 1e4` has to be resolved — so this is by far the
cheapest question in the campaign.


## RESULT (2026-07-27)

| `N` | 216 | 288 | 360 | 512 | 720 | 1440 |
|---|---:|---:|---:|---:|---:|---:|
| `t½` | 116.31 | **102.36** | 94.47 | 85.59 | 78.89 | 60.29 |

**Yang's `t½ = 107.42` falls between `N=216` and `N=288`** — it equals PARTIES at an effective
`N = 257`, against Yang's stated **288²** base grid for velocity and temperature (10.8 % apart).
Their 5× refinement to 1440² covers only salinity and the phase field, and the `Sm=0` case has no
salinity, so 288² is what sets their heat transport. PARTIES at 288² gives 102.36 vs 107.42 —
**4.7 %**, against ≈2.8 % digitization uncertainty.

The ladder is predictive **within the coarse range**: it forecast `N720` at 79.5 (got 78.89) and
`N288` at 101–102 (got 102.36). It is *not* predictive beyond it — see the extrapolation warning
below.

**Provenance PASSED (2026-07-27).** The `1440²` point is real. The re-run in
`../Yang_fresh_provenance/` (current binary, 128 ranks vs the original 256) reproduces the
completed run to **1.7e-12 at t=2**, 1.1e-5 at t=4, 5.8e-5 at t=6, against a 7.2e-3 discriminating
size at t=6 — a 125x margin. Binary and rank-count independence both confirmed, which also clears
the lineage under the salty production run.

⚠ **`t½(∞) ≈ 63–64` is NOT a converged limit — do not quote it.** Successive slopes
`dt½/d(1/N)` are 12058 / 11369 / 10764 / 11879 across the coarse rungs, then **26800** for
720→1440. The fit was calibrated on 216–720 and extrapolated a factor of two beyond it. The
ladder is not converged anywhere in the sampled range, 1440² included. The slope jump may be the
sidewall layer going unsteady between N=720 and N=1440 — hypothesis, not finding. PARTIES'
converged answer is `≤ 60.29`, i.e. further from Yang, which strengthens the headline.

⚠ PARTIES at `N≈258` has a band ≈23 % of `δ_T`; Yang keeps a sharp phase field and coarsens only
`u`,`T`. Same `t½`, possibly different mechanism — roadmap A.5.5 variant 3 must separate them.

Runs live in `/anvil/scratch/x-mjalabert/YangFresh_07272026/`.
Analyse with `python3 analyze_ladder.py`.

## Design

`make_inputs.py` generates `N360/`, `N512/`, `N720/`. Every physical input is byte-identical to
the completed `Sm=0` run (`Yang_production/parties.inp` with `cbd2 = cbd5 = 0`). Only the grid and
the four numbers tied to it change:

| `N` | `Δx = 1/N` | `Cn = melt_band_eps = 0.75Δx` | `Pe_CH = 0.9/Cn` | `dt₀` | ranks | est. cost to `t=140` |
|---:|---:|---:|---:|---:|---:|---:|
| 360 | 2.7777778e-3 | 2.0833333e-3 | 432.0 | 4.0e-4 | 36 | ≈ 110 core-h |
| 512 | 1.9531250e-3 | 1.46484375e-3 | 614.4 | 2.8e-4 | 64 | ≈ 315 core-h |
| 720 | 1.3888889e-3 | 1.0416667e-3 | 864.0 | 2.0e-4 | 144 | ≈ 875 core-h |
| 1440 | 6.9444444e-4 | 5.2083333e-4 | 1728.0 | 1.0e-4 | 256 | **already run**, `t½ = 60.287` |

Because `Cn = 0.75Δx`, refining the mesh also sharpens the diffuse band — this is a
diffuse-interface convergence study, not merely a mesh study.

`NConc = 2` with `s ≡ 0` is kept (not `NConc = 1`) so the ladder is exactly comparable to the
completed `1440²` point, which carried the inert salt field.

`t → 140` brackets both half-melt times (PARTIES 60.3, Yang 107.4).

## Run

Copy the deck (`parties.inp`, `jobscript.sh`, **`stop.inp`**) and the binary to a scratch
directory, then submit:

```bash
RUN=/anvil/scratch/$USER/YangFresh/N512
mkdir -p $RUN && cp N512/{parties.inp,jobscript.sh,stop.inp} $RUN/
cp $HOME/PARTIES/PARTIES/parties $RUN/
cd $RUN && sbatch jobscript.sh      # grids are independent; submit in any order
```

The binary must be built with the Stage-A flags (`EOS_NONLINEAR`, `CONC_VOF_PHASEWEIGHTED`,
`PHASE_CHANGE`, `ICE_PENALIZATION`, `TWOD_CARTESIAN`, `PERIODIC_Z_NOSLIP_BOX`).

**Binary provenance note.** The current binary (sha `6a8fa6eb…`, built 2026-07-21 03:16) is *not*
the one that produced the completed `1440²` `Sm=0` point (sha `aac878bf…`). It was rebuilt after
the Le=100 salt-operator work added the `yang_salt_transport` path. That path is a runtime flag
defaulting to **0**, is not set in these decks, and every use site in `Conc.c` is guarded by it,
so the physics path is the legacy one. Confirm empirically rather than assume: the early diffusive
phase is grid-independent physics, so all four grids must collapse onto the same 1-D Stefan
similarity curve before convection sets in. If they do not, the binary changed something.

### Two deployment gotchas — both cost a failed job on 2026-07-27

1. **`stop.inp` is mandatory.** PARTIES reads it every iteration and calls
   `Display_assert_error` if it is missing (`Temporal_int.c:495`), aborting at iteration 2. It is
   now generated with each deck; do not drop it when copying.
2. **Module names.** Use `gcc/11.2.0 openmpi/4.1.6 fftw/3.3.8` and launch with `srun`. The
   `gnu14/openmpi5/hdf5` names inherited from `Yang_production/jobscript.sh` do **not** exist on
   Anvil, and with `set -e` the job dies at the module load before the solver ever starts. HDF5 and
   HYPRE resolve from `~/Software/PARTIES_Libs` through the binary's RPATH and need no module —
   verify with `ldd parties`.

Analyse with `Yang_production/analyze_yang_production.py` and compare `t½` against
`Yang_production/yang_fig3a_reference.csv` (`Sm_g_per_kg = 0`).

## Precommitted decision rule

Let the ladder be `t½(360), t½(512), t½(720)` against the existing `t½(1440) = 60.287`.

- **Spread < 5 %** ⇒ converged at 360² already. Resolution is *not* the cause; the 1.78× is a
  model or reference discrepancy. Go to roadmap §A.5.5, and use this ladder as the cheap sandbox
  for every model variant.
- **`t½` rising with `N` toward ≈107** ⇒ the `1440²` run is under-resolved. Richardson-extrapolate
  before deciding whether a `2880²` freshwater confirmation (≈ 7 000 core-h) is worth it.
- **`t½` falling with `N`** ⇒ refinement moves away from Yang. Strong result; go to §A.5.5.

Record at matched `V/V₀` on each grid: interface superheat, liquid-side heat flux, `max_speed`,
and **mean interface length** (a scalloped front has more area than a flat one — the most likely
route by which resolution changes the melt rate).

## Context

`max_speed` in the completed `1440²` `Sm=0` run is 0.122–0.136 free-fall units over `t ∈ [4,60]`;
Yang's figure 4(b) velocity profiles peak at ≈0.085. The velocity scale is therefore **not** off
by 1.78 — the discrepancy is in heat delivery per unit velocity. See §6 of
`YANG_REFERENCE_CORRECTION_REPORT.md`, which also notes that a textbook sidewall-Nusselt estimate
(`t½ ≈ 68–77`) sits nearer PARTIES than Yang, so it is not safe to assume PARTIES is in error.
