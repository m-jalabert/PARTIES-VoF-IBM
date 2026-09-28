# A.4 Yang production run — paper-grounded expected results

> **TARGETS CORRECTED 2026-07-27 — §1.1 below was wrong.** Yang figure 3(a)/4(a), from which the
> `t½ ≈ 118` and `V(200)/V₀ ≈ 0.32` targets were digitized, is the **`ΔSv = 0`** series, not
> `ΔSv = 5`. This case is `Sm = 5, ΔSv = 5`, whose only published datum is the panel-(b) marker
> `f̄/f̄₀ = 0.497` ⇒ `t½ ≈ 216`. The error is visible in the two rows of §1.1 below that mix the
> panel-(a) `Sm`-curve ordering with the panel-(b) `ΔSv = 5` series. Corrected targets are in
> §1.1a; provenance and three independent proofs are in
> [YANG_REFERENCE_CORRECTION_REPORT.md](../../../../../docs/md_files/NASA-JPL/YANG_REFERENCE_CORRECTION_REPORT.md).
>
> **Post-run status (2026-07-12, re-scored 2026-07-27): overall FAIL, but in the opposite
> direction.** Morphology passed. `t½ = 179.848` against the corrected target `≈216` is **1.20×
> too fast**, not 52 % too slow. The companion `Sm=0` run (`t½ = 60.287` against the directly
> digitized `107.42`) is **1.78× too fast** and is now the campaign's largest discrepancy. See the
> [production analysis](../../../../../docs/md_files/NASA-JPL/YANG_PRODUCTION_ANALYSIS_REPORT.md).

**Case:** 2-D Cartesian reproduction of the reference case of Yang, Howland, Liu, Verzicco &
Lohse (2023), *Ice melting in salty water* (JFM 969, R2,
doi:10.1017/jfm.2023.582; arXiv:2302.02357; `docs/md_files/Yang2023.pdf`):
`H = 5 cm`, `ΔT = 20 K`, `Sm = 5 g/kg`, `ΔSv = 5 g/kg`, `RaT = 10⁷`, `Pr = 10`, `Sc = 1000`
(`Le = 100`), `St = L/(cpΔT) = 4`. Grid `1440 × 1440 × 1` in `TWOD_CARTESIAN`
storage-slab mode. Ice slab of
thickness `0.1H` on the east wall; quiescent start; run through the paper's comparison window,
`t/t_ff = 200`.

**Reference hierarchy:** numerical values below are read from the published JFM paper, which is
the ground truth. Published figure numbers are used first; the arXiv letter has the volume curve
as Fig. 3 and the layer plot as Fig. 2. The PARTIES run should be compared with the paper's
**2-D** result; `z` is a one-cell periodic storage slab, not a physical third dimension.

**Known method differences (JFM-verified 2026-07-12):** Yang's Figure 1 grid study is for this
exact `Sm = 5`, `ΔSv = 5` case. It chooses base `nx = 288`, with salinity and phase field fixed
at five-times refinement, i.e. **288² for velocity/temperature and 1440² for salinity/phase**.
The paper's 432²/2880² example applies specifically to `RaT = 10⁸`, `RaS = 2×10¹⁰`, not this
reference case; its separate 3-D example uses 288×288×144 and 864×864×432. Their phase field is
an Allen–Cahn equation with `ε = grid spacing`, `D = 1.2κ`, `εΔT/γ = 1`; velocity is damped with
penalty `η = dt` plus **direct forcing to zero for solid fraction > 0.9**. PARTIES instead uses
one 1440² grid (matching their refined S/phase resolution and 5× finer than their u/T base grid),
its validated CH
parameters (`Cn=0.75Δx`, `Pe_CH=0.9/Cn`), and a fixed implicit Brinkman time `darcy_tau=1e-4`
without the direct-forcing cutoff (≈4× stronger drag at the band edge for a typical dt). The
liquidus slope `m = 0.056 °C/(g/kg)` is stated explicitly in the JFM text, confirming
`liquidus_slope = 0.014`. These are controlled implementation differences, not values claimed
to be identical to the paper.

**Unit conversion:** the free-fall velocity is `U_ff = √(g·(CbΔT²/2ρ₀)·H) ≈ 3.3 cm/s`, so one
free-fall time is `t_ff = H/U_ff ≈ 1.5 s`; `t = 200` corresponds to ≈ 5 minutes of physical
melting.

---

## 1. Primary benchmark criteria (the Stage-A validation)

These two comparisons are the A.4 pass criteria (roadmap §A.4). The A.3–A.3d gates isolate the
individual implementations, but they do not prove that PARTIES' CH mobility/melt coupling and
Brinkman treatment are quantitatively identical to AFiD in the full problem. A discrepancy must
therefore be diagnosed before assigning it solely to resolution or geometry.

### 1.1 Ice volume V(t)/V₀ vs published Yang Fig. 4(a) (arXiv Fig. 3a)

`V(t) = ∫(1−F) dV`; compute `V₀` from the stored `t=0` field (nominal sharp-slab value
`0.1·Ly·Lz = 0.1·Lz`, where `Lz` is the bookkeeping slab thickness). Read from the published figure
(ΔSv = 5 g/kg family; the Sm = 5 curve is the light-blue one):

| Quantity | Expected value | Notes |
|---|---|---|
| Shape of V(t)/V₀ | smooth and monotonic, but **do not assume a constant slope** | digitization gives slope −0.003827 over t=40–100 and −0.002292 over t=120–200 |
| Half-melt time t½ (V = 0.5·V₀) | **≈ 118 t_ff** (figure digitization uncertainty ≈ ±3) | the Sm = 5 curve crosses 0.5 last among Sm = {0, 5, 10, 15} |
| V(200)/V₀ | **≈ 0.32** (roughly 0.30–0.34 by line thickness) | Sm = 5 is the slowest-melting of the four plotted curves |
| Mean melt rate f̄ = 1/t½ | **≈ 0.0085 t_ff⁻¹**; `f̄/f̄₀ ≈ 0.5` | Sm = 5 is close to the ΔSv = 5 minimum near Sm ≈ 6 in published Fig. 4(b) |

Project acceptance: `t½` within 10 % and `V(200)/V₀` within about 0.05 of the digitized 2-D
curve. These tolerances are project choices, not uncertainty bounds stated by Yang et al. The
ordering is a sharper later check: another `Sm` point must reproduce the non-monotonic trend.

### 1.1a CORRECTED targets (2026-07-27) — supersede the table above

The table in §1.1 reads the panel-(a) curve family (`ΔSv = 0`, four `Sm` values) as if it were the
panel-(b) `ΔSv = 5` series. It is not. Digitized from the archived arXiv vector figure by
`digitize_yang_fig3.py` (data: `yang_fig3a_reference.csv`, `yang_fig3b_reference.csv`):

| Case | Yang `t½` | Yang `V(200)/V₀` | source |
|---|---:|---:|---|
| `Sm=0` (`ΔSv=0` necessarily) | **107.42** | 0.2483 | fig 3(a), measured directly |
| `Sm=5, ΔSv=0` | 118.42 | 0.3142 | fig 3(a) |
| `Sm=10, ΔSv=0` | 112.21 | 0.2837 | fig 3(a) |
| `Sm=15, ΔSv=0` | 101.84 | 0.2245 | fig 3(a) |
| **`Sm=5, ΔSv=5` — THIS CASE** | **≈216** | not published | fig 3(b) rank 3: `f̄/f̄₀ = 0.497` |

`f̄ = 1/t½` is confirmed verbatim from the paper source, so `f̄₀ = 1/107.42`.

**Whenever a target is quoted from now on, state `ΔSv` with it.** The whole error was possible
because the campaign wrote "the Sm = 5 curve" without its `ΔSv`, and the paper plots two different
`Sm = 5` cases in the two panels of the same figure.

`yang_fig4_sm5_reference.csv` is **superseded**: it holds the `Sm=5, ΔSv=0` curve under a
`Sm5_DeltaSv5` label. Retained for provenance only; do not score against it.

### 1.2 Melt-front layering vs Huppert & Turner (published Yang Fig. 3b; arXiv Fig. 2b)

With the EOS of the paper, at `T∞ = 20 °C`, `S∞ = Sm = 5 g/kg`:

- horizontal density contrast: `Δρ_T = ρ(0, S∞) − ρ(T∞, S∞) = 1.595 kg/m³`,
- ambient stratification from the full EOS at `T∞=20 °C`, `S∞=5 g/kg`:
  `|dρ/dz| ≈ 72.26 kg/m⁴` (i.e. `ρ₀⁻¹|dρ/dz| ≈ 7.23·10⁻⁴ cm⁻¹` for `ρ₀≈1000 kg/m³`),
- Huppert–Turner prediction: `h = 0.65·Δρ_T/|dρ/dz|` ⇒ **h ≈ 1.43 cm = 0.287H**, i.e.
  **3–4 layers** stacked over the domain height.

Measure h as the mean vertical spacing of the scallops in the melt-front contour
(`F = 0.5` isoline) over the middle of the domain, once layers are developed (t ≳ 50).
Acceptance: the point (`7.23·10⁻⁴ cm⁻¹`, `h·ρ₀/Δρ_T ≈ 900 cm`) must sit inside the
experimental/numerical scatter of published Fig. 3(b); the H&T coefficient is 0.65 ± 0.06, so ±10 %
plus half a layer of measurement ambiguity is expected.

---

## 2. Paper-supported qualitative evolution

The paper supports the sequence below but does not publish precise onset times for the reference
case. Any times used during monitoring are project checkpoints, not ground-truth thresholds.

1. **Early stage:** near-diffusive melting; a thin fresh, cold meltwater film forms on the front
   and rises (salinity-driven buoyancy dominates in the film because ΛT = 1.75 > 1).
2. **Layer-development stage:** outside the film, water cooled at the front (heat diffuses 100× faster
   than salt) sinks until it reaches neutral buoyancy in the ambient stratification —
   vertically stacked convection rolls appear.
3. **Developed stage:** the rolls sculpt a **layered/scalloped melt front** with ~3–4 layers
   (criterion 1.2). Cold fresh water accumulates under the lid, so the very top of the ice
   melts *slowest*: expect a **local maximum of ice thickness at the top boundary** — this
   is the signature feature of the ΔSv = 5 reference case.
4. The reference case sits at `ΛT = 1.75`, `ΛS = 1` — near the regime II/III boundary of
   Yang Fig. 4(c), which is why its melt rate is close to the minimum of the ΔSv = 5 curve.
5. The paper shows near-front velocities of order `0.1 U_ff` for selected cases in published
   Figs. 5(b)/6, but it does not provide a numerical velocity target for this exact
   (`Sm=5`, `ΔSv=5`) reference case. Use non-zero layered circulation as a qualitative check,
   not `0.05–0.1` as an acceptance band.

---

## 3. Run-health diagnostics (should hold the entire run)

All physical walls are no-flux and the storage slab is periodic, so the box is closed — the same total-conservation
statements certified in the gates apply here **with the flow on**:

| Check | Expectation | Certified reference |
|---|---|---|
| Salt budget `∫s dV` | constant to ≤ 10⁻⁸ relative (gates achieved ~10⁻¹³–10⁻¹⁶) | A.3b, A.2 |
| Enthalpy `∫(θ + F/stefan) dV = ∫(θ + 4F)dV` | constant to ≤ 10⁻⁶ relative over diagnostic windows (`stefan=cpΔT/L=0.25`; paper's `St=L/(cpΔT)=4`) | A.3, A.2 |
| Salt in deep ice (x > front + 25Δx) | ≤ 10⁻⁶ | A.3b, A.2 |
| F bounds | F ∈ [0, 1] exactly (conservative clip active) | A.3 |
| Deep-ice speed | ≪ free-fall unit; penalized to O(10⁻⁴) or below | A.2, A.3d |
| dt (adaptive) | `≤ 5·10⁻⁴`; actual value is a run diagnostic, not specified by the paper | input CFL/cap |
| Divergence residual | < 10⁻⁶ each projection | A.2 |
| max|u| | not quiescent after spin-up; exactly-zero velocity signals inactive buoyancy | A.3c and paper flow fields |

The explicit melt-source stability bound `dt ≲ Pe_T·ε·4Cn ≈ 10⁻²` is far above the CFL
limit here and cannot bind.

---

## 4. What a discrepancy would mean

| Symptom | Likely cause (gates exclude the rest) |
|---|---|
| t½ off by ≫ 10 % but layering correct | inspect CH mobility/melt mapping and Brinkman rigidity; `mu2=1` intentionally matches the paper's constant viscosity |
| No layers / wrong layer count with correct Δρ_T, dρ/dz | check salt-boundary-layer resolution and `max|s|` near the front; 1440 matches the refined S/phase grid of the exact paper case, but PARTIES still needs its own convergence check |
| Salt budget drifting | new integration-level bug (all masked-operator paths were machine-conservative in A.2/A.3b) |
| Melt front advancing uniformly with no top-thickness maximum | stratification not active — check init 31 profile (`cbd2 = 0.5` top, `cbd5 = 1.5` bottom) and `eos_betaS` |
| V(t) stalls relative to the reference curve | compare against published Fig. 4(a); do not label it expected without that comparison |

---

## 5. Practical notes

- **Outputs:** `output_time_interval = 2.0` gives about 101 output times through `t=200`.
  The fixed 2-D storage-slab run writes roughly one quarter of the failed four-cell slab output:
  about **0.35 GB per `Data_*` plus 0.25 GB per `Resume_*`**, or about **60 GB total**.
  `Resume.h5` is only a symlink. Run from
  `/bigscratch`, not home, and confirm quota first.
- **Wall time:** at most `dt=5e-4` means at least 400,000 steps to `t=200`; smaller CFL steps
  can raise this above one million. This will probably span several SLURM
  windows: set `resume = 1` in `parties.inp` after the first job and resubmit (restarts
  from the latest `Resume.h5`).
- The production input stops at **t = 200**, the end of the published volume-curve window.
  Extending beyond 200 is optional and is not part of this benchmark acceptance.
- `stop.inp` set to 1 stops the run gracefully at the next step.
- Half-melt happens near `t≈118`. Inspecting morphology around `t=50` is a useful project
  checkpoint, but the paper does not identify `t=50` as the onset of layering.
