# ECCO_TESTS — test cases of the Yang benchmark → ECCO sediment-from-ice campaign

Reorganized by stage on 2026-09-28 (`MOVES_20260928.tsv` lists every old → new
path; a tarball of the previous layout is at
`/anvil/scratch/x-mjalabert/ECCO_StageB3/backups/ECCO_TESTS_before_reorg_20260928.tar.gz`).
The controlling plan is the
[implementation roadmap](../../../docs/md_files/NASA-JPL/YANG_ECCO_IMPLEMENTATION_ROADMAP.md).
Raw run output lives on scratch, never here.

| Folder | Roadmap stage | Contents |
|---|---|---|
| `StageA_Yang/` | Stage A (closed 2026-07-30) | Yang production and provenance runs, fresh/salty resolution ladders, clip audit, Le=100 salt-operator gate (+ its summary CSV/MD), Melting-RB discriminator |
| `StageB1_development/` | B.0 / B.1 | Non-invasiveness checks, Stage-A no-impact gate, release and tracer development, B.1.4 sediment-clip arms, polar EOS fits |
| `StageB2_closure/` | B.2 (closed 2026-09-08) | `B2_closure/` (closure campaign, validated header `Boundary.validation.h`, `build_validation.py`, `analyze_profiles.py`), `Settle_validation/`, `StageB_presentation/` (figures/videos) |
| `StageB3_scenarios/` | B.3 | `P1a_feasibility/` (first P1, job 20906327: no release, superseded) and `P1b_shelf/` (released fine-sand grain beneath an Antarctic ice-shelf base, melting ON) |

`StageB_release2d`, `Yang_clip_audit` and `Yang_salty_provenance` are symlinks to
`/anvil/scratch/x-mjalabert/ECCO_StageB_output/`. The roadmap's historical
`ECCO_TESTS/Yang_checkout/` reference predates this folder and has no copy here.

Scripts locate each other relative to their own folder, so keep the stage
folders as siblings. `build_validation.py` finds the solver root by walking up
to the nearest folder that holds `src/` and the `Makefile`.
