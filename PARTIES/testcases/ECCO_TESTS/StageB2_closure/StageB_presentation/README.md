# Stage B — presentation figures and videos

Twelve figures and two videos summarising **STAGE B — ECCO sediment-from-ice study**,
for the talk that precedes
[§B.3 Full ECCO production inputs](../../../../../docs/md_files/NASA-JPL/YANG_ECCO_IMPLEMENTATION_ROADMAP.md).

The controlling numerical record is still
[the B.2 closure report](../../../../../docs/md_files/NASA-JPL/YANG_ECCO_B2_CLOSURE_REPORT.md)
and the Stage-B sections of the roadmap. **Nothing here is new evidence**, with one
explicitly-marked exception (figure 11d, below). These scripts re-derive the recorded
numbers from the run directories so the slides can be regenerated and checked.

---

## Contents

```
stageb_common.py        paths, run registry, loaders, palette, plot style
make_all.py             regenerate everything
fig01..fig12_*.py       one script per figure
video1_3d_lifecycle.py  3-D movie
video2_slices.py        four-field mid-plane movie
figures/                *.png (200 dpi) and *.pdf (vector) — committed outputs
videos/                 *.mp4 (h264, 12 fps, ~11 s each)
data/                   small derived artefacts the scripts cache or emit
```

### The figures

| # | File | What it is for |
|---|---|---|
| 1 | `fig01_stageb_map` | Opening slide: the four Stage-B blocks, the six B.2 gates as measured, and the constraint B.2 hands to B.3 |
| 2 | `fig02_configuration` | The box to scale, the initial state read back from `Data_0.h5`, the whole deck read back from `parties.inp`, and which of those numbers are synthetic |
| 3 | `fig03_eos_polar` | B.1.1 / B.1.1a: the polar EOS re-fit, the S ≈ 24 crossover, haline-over-thermal dominance, and why the linear EOS is the safer choice |
| 4 | `fig04_release_mechanism` | B.1.3 hold-in-ice release, the `F_release` calibration, and the uncompensated `C_L` clip that the meltwater tracer exposed (B.1.4) |
| 5 | `fig05_lifecycle` | The headline result: locked → released → settling → landed, with the four motion gates |
| 6 | `fig06_conservation` | Heat / salt / tracer budgets, ice rigidity, and the negative controls |
| 7 | `fig07_reproducibility` | Where bit-identity is claimed and where only roundoff agreement is, plus the non-invasiveness contract |
| 8 | `fig08_conditioning` | The one quantity that does **not** converge — the escape from the ice — and the measurements that explain it |
| 9 | `fig09_resolution` | Settling refinement, the Darcy isolation controls, and an explicit list of what is *not* claimed |
| 10 | `fig10_fields` | Mid-plane montage: temperature and meltwater anomaly at four times |
| 11 | `fig11_mixing` | Where the meltwater goes (in-solver profiles), plus **one flag for B.3** |
| 12 | `fig12_cost` | What Stage B cost, and what B.3 needs against what is left |

A suggested talk order is 1 → 2 → (video 1) → 5 → 4 → 6 → 10 → (video 2) → 11 → 7 → 8 → 9 → 3 → 12.
Figures 8, 9 and 11d are the ones a reviewer will ask about; they are deliberately blunt.

### The videos

Both are made from the accepted B.2 coupled gate, **job 20492853**
(`B2closure_coupled16_rampfix_20492853`), 128 frames at 12 fps.

| File | What it shows |
|---|---|
| `videos/video1_3d_lifecycle.mp4` | Left: marching-cubes isosurface of the **actual ice** with the grain; right: meltwater anomaly on the plane through the grain; below: a timeline of the grain height |
| `videos/video2_slices.mp4` | Temperature, meltwater anomaly, vertical velocity and speed on that same plane |

**Time convention, stated on every frame.** PARTIES wrote 33 field outputs, one every
2.0 time units. The Eulerian fields are linearly blended between consecutive outputs so
the movie is watchable; the grain's position and velocity come from `mobile.dat` at full
step resolution, so the grain moves exactly as it did in the run. During the fall the
grain covers ≈ 0.65 diameters per saved output, so the blended *fields* around it are
interpolated, not resolved — read the fall as a trajectory, not as a resolved wake.

**What is drawn as ice.** Everywhere in this folder "ice" means the project's own
`real_ice`: `max(0, 1 − C_L − C_S)` masked to `C_S < 0.05`. Inside the IBM's excluded
support `C_L` is an extension value, and `1 − C_L − C_S` there is not ice. Scalars inside
that support are masked out of every plot for the same reason.

---

## Regenerating

On Anvil, from this directory, with the anaconda python (NumPy, h5py, Matplotlib, SciPy,
scikit-image, and `gsw` for figure 3 only):

```bash
python make_all.py                 # the twelve figures, ~45 s total
python make_all.py --videos        # and both movies, ~4 min more
python make_all.py --only fig08    # just one
```

Videos need an ffmpeg. `stageb_common.use_ffmpeg()` takes `$FFMPEG_BINARY`, else the
binary bundled with `imageio-ffmpeg`, else whatever is on `PATH` — so
`module load ffmpeg/4.2.2` also works. No MPI and no SU are needed: everything here reads
saved output.

Figure 11d scans two runs' full field output (66 × 95 MB) and caches the result in
`data/low_ice_scan.json`; pass `--rescan` to redo it.

---

## Where the numbers come from

Three sources, and the scripts keep them separate on purpose:

1. **Run directories on scratch** — `/anvil/scratch/x-mjalabert/ECCO_StageB/…`, registered
   by name in `stageb_common.RUNS`. `mobile.dat`, `release.dat`, `ecco_profiles.csv`,
   `Data_*.h5`, `Particle_*.h5`, `parties.inp`.
2. **Recorded audits** in [`../B2_closure`](../B2_closure/README.md) — the JSON files that
   the closure campaign produced, plus its analysers. Figures 6 and 8 *import*
   `analyze_profiles.py` and `refine_timesteps.py` rather than re-implementing their
   budgets and metrics, so a figure cannot silently disagree with the audit.
3. **Prose-only values** from the roadmap and the closure report, collected in
   `stageb_common.ROADMAP_NUMBERS` with the section each came from. These are the values
   whose underlying runs are not on scratch any more (B.0 part 2, the 2-D B.1.4
   discriminator, the `F_release` calibration, the B.1.1a extrapolation table, the
   allocation and forward-cost tables).

Every figure prints its sources in its own caption strip, and its docstring says which of
the three it used.

---

## Two things to read carefully before presenting

**1. Figure 8 is not a failure slide, and should not be presented as one.**
The matched-ramp timestep comparison fails its debugging tolerances, and that is the
correct outcome: the grain's escape from the ice is a near-cancellation (the IBM reaction
is 97.9 %, 99.0 %, 99.9 % of the buoyant weight as `dt` halves), so its timing drifts
instead of converging. Chaos, Darcy leakage and a dt-dependent melt rate were each tested
and refuted before conditioning was accepted. The settling dynamics themselves are
converged to 0.08 % at `max_dt = 0.01`; the hard ice-penalization threshold costs about
8× of that, and disabling it is a trade-off, not a fix.

**2. Figure 11d is an observation made while preparing these slides, not a recorded gate.**
Scanning the saved fields of the two twins shows that the deck with
`VOF_DIFFUSE_ICE_PENAL_THRESHOLD` **active** grows downward ice tongues after release —
peaking at 0.33 volume units, about 0.9 % of the initial ice, at t ≈ 46, then melting away
by t ≈ 62 — while the threshold-disabled twin, run with the same deck to the same time,
has **exactly zero** actual ice below y = 3 at every saved output. The measurement uses
the project's own `real_ice` definition, so it is not the IBM support being miscounted,
and a hard 0/1 Darcy flip at ice fraction 0.9 is a plausible mechanism (a damped cell
stops advecting heat, stays cold, and keeps freezing). It has **not** been through the
campaign's audit chain: no controlled twin at matched melt history, no resolution or
timestep check, no confirmation that it is not an artefact of the `> 0.9` counting
threshold. Present it as a question for B.3 — *does the hard threshold need to stay on?* —
not as a result.

Related, and already recorded: the closure report notes that the same threshold flips are
discrete O(dt) events in the grain's wake, degrade settling convergence about 8×, and that
disabling it delays release by 10.3 time units and worsens the ice-rigidity gate.

---

## Caveats carried from the campaign

- Everything in B.2 uses **synthetic debugging parameters**. The salt field is passive
  (`eos_betaS = 0`) so it exercises conservation without pretending to be the polar EOS;
  `Pe`, `St`, `Cn`, `τ`, `max_dt` and `d/Δx = 16` were chosen to make the shake-out
  affordable. None of them is a dimensional production value.
- The ΔT and S_m in figure 3 are the two **placeholder** ECCO scenarios. The dimensional
  TEOS-10 fit is scenario-independent and final; the nondimensional block must be
  regenerated once the physical case is fixed.
- Figure 9 says out loud what the settling runs do *not* establish: no terminal velocity,
  no Schiller–Naumann validation, no confinement correction.
- No mixing efficiency is inferred anywhere. An open melting box needs its source and
  boundary BPE accounting, and the run's actual EOS, before η means anything.
- The three settling control jobs in figure 9c are recorded **TIMEOUT**, not completed
  runs; their shared output through t ≈ 11 is what supports the comparison.
