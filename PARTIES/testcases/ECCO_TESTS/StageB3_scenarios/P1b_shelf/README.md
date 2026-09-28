# P1b — released fine-sand grain beneath an Antarctic ice-shelf base

The full design, pre-flight criteria and run record are in
[YANG_ECCO_B3_P1b_SHELF.md](../../../../../docs/md_files/NASA-JPL/YANG_ECCO_B3_P1b_SHELF.md).

| File | Purpose |
|---|---|
| `prepare.py` | Regenerates everything below from the physical record (needs `gsw`; use the Anvil anaconda 2025.06 python) |
| `case.json` | Three-equation state, EOS/liquidus fits, scaling, geometry, cost, gate reference |
| `parties.inp`, `p_mobile.inp`, `p_fixed.inp`, `stop.inp`, `Boundary.scenario.h` | P1b production deck |
| `ctrl/` | Matched ice-only control (no grain) |
| `smoke/` | Full-size deck to t=0.2 (initial-condition and cost check) |
| `gate/` | 1-D kernel-gate variants and `reference.json` (exact binary-Stefan solution) |
| `job_gate_smoke.sh`, `analyze_gate.py`, `check_smoke.py` | Pre-flight job and its checks |
| `job.sh` | Production job (`CASE=ctrl` for the control) |

Stage a copy under `/anvil/scratch/x-mjalabert/ECCO_StageB3/P1_Attempts/P1b_shelf/`
and submit from there; each job writes a job-unique run directory.
