# Stage-A Intermediate Validation Report — 1D Binary (Salty) Stefan Gate (A.3b)

**Date:** 2026-07-03
**Scope:** The A.3b intermediate physics gate of the
[Yang–ECCO roadmap](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) — fresh ice melting into warm **salty**
water against the exact coupled similarity solution (two-phase Neumann + meltwater-dilution
boundary condition + liquidus depression) — run to completion in
`PARTIES/testcases/1DStefanSalty/`. This is the companion to the
[A.3 Stefan gate report](STAGE_A_STEFAN_GATE_REPORT.md) and the **last physics test before the
A.4 Yang validation**.
**Verdict:** **PASS on all eight roadmap criteria, on the first attempt, with zero code
changes** — the melting kernel built and debugged for A.3 generalized correctly to the coupled
salt problem.

---

## 1. Headline results

1D binary Stefan gate (768 × 4 × 4 thin strip, `Lx = 1.5`, `dx = 1/512`, interface at
`x = 0.75`, `Pe_T = 100`, `Pe_S = 1000` (Le = 10), `stefan = 0.5`, `liquidus_slope = 0.5`,
`θ_∞ = 1`, `θ_ice = 0`, `s_∞ = 1`, all-no-flux box, 4 MPI ranks, `dt = 2·10⁻⁴` constant,
t → 4.5, 22 500 steps, ≈ 10 min wall time):

| Roadmap criterion (§A.3b) | Requirement | Result |
|---|---|---|
| 1. Front law `X(t) = 0.75 + 2λ√(t/Pe_T)` | λ = 0.25554 ± 3 % | **λ = 0.25320 (0.92 % error)**, fit RMS 4.0·10⁻⁵ ≈ Δx/49 |
| 2. Temperature similarity (two-sided) | max error ≤ 0.05 outside the band | 1.29–1.66·10⁻² at t = 1.8, 2.7, 3.6, 4.5 (**decreasing** in time) |
| 3. Salinity similarity (dilution profile) | max error ≤ 0.05·s_∞ | 3.7–4.5·10⁻³ (decreasing in time) |
| 4a. Interface salinity `s(X)` | s_i = 0.1722 ± 15 % | mean 0.1646 (4.4 %), converging 0.1596 → 0.1666 |
| 4b. Interface temperature `θ(X)` | θ_i = −0.0861 ± 0.03 | mean −0.0747 (0.011), converging −0.0666 → −0.0780 |
| 5. Salt conservation `∫s dV` | drift ≤ 10⁻⁸ relative | **2.9·10⁻¹³** (0.75 exactly, all outputs) |
| 6. No salt in ice (`x > X + 25Δx`) | ≤ 10⁻⁸ | ≤ 3.6·10⁻¹¹ (initial tanh tail; decays to 10⁻²²⁴) |
| 7. Enthalpy conservation `∫(θ + F/St) dV` | drift ≤ 10⁻⁸ relative | **3.0·10⁻¹⁰** (2.25, all outputs) |
| 8. Quiescence | max\|u\| < 10⁻¹⁰ | **max\|u\|,\|v\|,\|w\| = 0.0 exactly**, entire run |
| MPI ([MPI] requirement) | rank-count independent | 1-rank vs 4-rank at t = 0.45: Δθ ≤ 2.1·10⁻¹⁵, Δs ≤ 1.0·10⁻¹⁵, ΔF ≤ 3.2·10⁻¹⁵, Δu = 0 |

The automated gate is `PARTIES/testcases/1DStefanSalty/analyze_salty_stefan.py` (exit 0 = pass);
it writes `stefan_salty_results.csv` (gate summary), `stefan_salty_timeseries.csv`,
`stefan_salty_profiles.csv`, and the `fig_*.png` figures (front law with broken-coupling
references, two-sided θ collapse, s dilution collapse, interface values, conservation,
raw profiles).

Because criteria 5 and 7 hold in an **all-no-flux box** (melting driven purely by the initial
liquid superheat), they are *total*-conservation statements — strictly stronger than the A.3
wall-flux closure, with no time-quadrature caveat.

---

## 2. The problem, its exact solution, and what is newly certified

The A.3 gate certified the melting kernel for pure water; a Yang discrepancy could still not have
been attributed between the salt couplings. This gate pins them down in isolation: liquid
(`θ = 1`, `s = 1`) for `x < 0.75`, fresh ice (`θ = 0`, `s = 0`) beyond, both phases effectively
semi-infinite over the run, `u ≡ 0`. Uniform initial states in each phase *are* the `t → 0` limit
of the similarity solution, so the front law holds essentially from the start (fitted virtual
origin `t₀ = +0.020`).

With `η = (x − 0.75)/(2√(αt))`, `α = 1/Pe_T`, `λ_s = λ√Le`, the three interface conditions close
the system for `(λ, s_i, θ_i)`:

```
liquidus:       θ_i = T_melt − Λ*·s_i
salt dilution:  (1/Pe_S)·∂s/∂x|_i = −Ẋ·s_i   ⇒   s_i = s_∞ / [1 + √π·λ_s·e^{λ_s²}·(1+erf(λ_s))]
Stefan:         λ√π e^{λ²}/stefan = (θ_∞ − θ_i)/(1 + erf(λ)) + (θ_ice − θ_i)/erfc(λ)
```

For the gate inputs: **λ = 0.25554, s_i = 0.17220, θ_i = −0.08610** (solved in
`analyze_salty_stefan.py`). Both Stefan terms *add*: the liquidus depression puts the interface
below both far fields, so salt makes ice melt **faster** (the road-salt effect) — the front-law
plot shows the measured points riding the exact curve, clearly above both broken-coupling
references.

**Discrimination.** The measured λ separates every plausible failure mode by far more than the
0.92 % measurement error:

| Coupling state | λ | Distance from measured 0.25320 |
|---|---|---|
| Everything correct | 0.25554 | 0.92 % |
| Solid-side conduction lost (`κ_r(F)` path broken) | 0.23391 | 7.6 % |
| Liquidus/salt coupling lost (`θ_L(s)` ignored) | 0.21688 | 14.3 % |
| Dilution BC lost (`s_i → s_∞`) | (θ_i → −0.5, front far too fast) | — |

**Code paths newly certified** (none of which A.3 exercised):

| Path | Where | How this gate checks it |
|---|---|---|
| Salinity-dependent liquidus `θ_L = T_melt − Λ*·s` | `VOF_DIFFUSE_compute_melt_rate` (`has_salt`, A.1.4) | λ and θ(X) → −0.086, not 0 |
| F-masked salt diffusivity (`kappa_ice_ratio = 0`) | `Conc_phase_kappa` both CN halves (A.1.3) | salt-in-ice ≤ 3.6·10⁻¹¹; ∫s exact |
| **Emergent meltwater-dilution BC** (G.5 claim) | no explicit BC — species conservation at the moving front | s(η_s) collapses on the exact dilution profile; s(X) → s_i |
| Solid-side conduction (`kappa_ice_ratio = 1`) | `κ_r(F)` in both CN halves (A.1.3) | two-sided θ collapse; the solid term is 12 % of the heat budget |
| `NConc = 2` field layout (0 = θ, 1 = s) + per-field `Pe`/`kappa_ice_ratio`/BC arrays | Yang configuration exactly | the whole gate |

Deliberately **not** exercised (no analytic quiescent-1D reference exists): `EOS_NONLINEAR`
buoyancy and `ICE_PENALIZATION` — covered by the A.2 shakeout diagnostics and the Yang comparison
itself.

---

## 3. Design choices that make the gate sharp

1. **Exaggerated liquidus** (`liquidus_slope = 0.5` vs Yang's 0.014): at the Yang value the
   coupling is a 1 % perturbation, indistinguishable from discretization error; at 0.5 it moves λ
   by 15 % and the gate becomes unambiguous. The melt law is linear in `θ − θ_L(s)`, so
   correctness at Λ* = 0.5 implies correctness at 0.014.
2. **Le = 10** (`Pe_S = 1000`): large enough that the salt layer is distinctly thinner than the
   thermal layer (√Le ≈ 3.2×) and the dilution structure is genuinely two-scale; small enough
   that the salt boundary layer (23 cells at the fit-window start, 69 at the end) stays well
   wider than the CH band (~4 cells). The Yang-scale Le = 100 steepening is a *resolution*
   question, already handled by the 1440² grid choice in A.4 — not a coupling question.
3. **All-no-flux box, superheat-driven**: no hot wall means no wall-flux quadrature in the
   budgets — salt and enthalpy conservation become machine-precision checks (10⁻¹³ / 10⁻¹⁰) that
   would expose any non-conservative term in the F-masked operators immediately.
4. **Both phases semi-infinite by construction**: `Lx = 1.5` with the interface at the middle
   keeps the liquid gap ≥ 4.1 and the (shrinking) solid gap ≥ 3.0 diffusion lengths through
   t = 4.5, so the far fields stay clean over the whole fit window (t ∈ [0.5, 4.5], front travel
   35 cells).
5. **`θ_ice = 0` with `θ_i < 0`**: the solid conducts heat *into* the depressed interface — the
   same sign the term has in the real Yang problem. This supersedes the "optional subcooled
   two-phase variant" of A.3 (the solid-side branch is now certified in its production-relevant
   direction).

---

## 4. Error accounting — why 0.92 % and why it is benign

All residuals behave as the **shrinking band-scale bias** expected of a diffuse-interface method,
not as a physics error:

- the origin-free instantaneous constant `λ_inst(t) = √(d(X−x₀)²/dt / 4α)` climbs monotonically
  0.2517 → 0.2542 across the fit window, approaching the exact 0.25554 from below as the
  boundary layers widen relative to the fixed band width (contrast A.3 iterations 2–3, where a
  *flat* λ_inst deficit diagnosed a genuine leak);
- the similarity errors fall in time (θ: 1.66 → 1.29·10⁻²; s: 4.5 → 3.7·10⁻³);
- `s(X)` and `θ(X)` converge monotonically toward `s_i`, `θ_i`; their offsets (4.4 %, 0.011) are
  the O(half-band × interface gradient) *sampling* uncertainty of reading a steep profile at the
  `F = 0.5` point — at t = 1 the salt gradient at the front is ≈ 4.4, so half a band (≈ 2Δx)
  subtends ≈ 0.017 in s, which is the whole observed offset;
- the residual λ deficit is the same O(`melt_band_eps`) superheat bias measured in A.3 (0.4 %
  there at ε = Cn ≈ 0.003 with Pe_T = 10), here compounded by the salt-in-band structure — both
  contract with `Cn`/`melt_band_eps`/grid refinement if a tighter number is ever needed.

No iteration was needed: the formulation-level defects (one-sided flux jump, |∇F| deposition,
non-conservative clip) had already been found and fixed by the A.3 gate, which was the point of
running the gates in this order.

---

## 5. MPI rank-count consistency

Same procedure as A.3: identical inputs re-run on 1 rank and compared field-by-field against the
4-rank production run at t = 0.45 (2 250 steps, `Data_3.h5`):

| Field | max abs. difference (1 vs 4 ranks) |
|---|---|
| θ (`Conc/0`) | 2.1·10⁻¹⁵ |
| s (`Conc/1`) | 1.0·10⁻¹⁵ |
| F (`VOF/C_L`) | 3.2·10⁻¹⁵ |
| u, v, w | 0.0 exactly |

The melt law is pointwise, the salt masking enters only through face-averaged coefficients in the
existing domain-decomposed CN operators, and the conservative-clip redistribution is
MPI-collective — so rank-count independence is by construction; this check confirms it.

---

## 6. Documentation updates applied

1. **Roadmap §A.3b added** (`YANG_ECCO_IMPLEMENTATION_ROADMAP.md`): full problem statement,
   exact solution, discrimination table, `Boundary.h` flags, inputs, eight pass criteria,
   resolution caveats — then annotated ✅ PASSED 2026-07-03 with the as-run numbers.
2. **Roadmap A.3 annotated** ✅ PASSED 2026-07-02 with headline numbers and a link to the A.3
   report; the stale "code prerequisite: A.2 shakeout passed" wording was corrected to the
   as-built order (A.3/A.3b ran on the reduced flag set; A.2 with the **full Yang flag set** is
   still pending and required before A.4).
3. **Scope/omissions blocks updated**: Stage A now keeps exactly **two** 1D physics gates; the
   A.4 stage-transition note and quick-reference tables mention A.3b (init 27 at
   `vof_slab_x0 = 0.5`; init 31 with `cbd2 = cbd5` ⇒ uniform salinity).
4. **A.3 report §6** updated: intermediate gate marked done, superseding the "optional
   subcooled two-phase variant".

---

## 7. Files touched

**No solver source files were modified.** The gate consumed only existing, already-validated
code paths.

| File | Change |
|---|---|
| `src/Include/Boundary.h` | comment-only: the phase-change flag block now documents that the current configuration serves both the A.3 and A.3b gates (flags themselves unchanged; `NConc` is a runtime input) |
| `testcases/1DStefanSalty/` | **new**: `parties.inp`, `p_mobile.inp`/`p_fixed.inp` (0 particles), `stop.inp`, `analyze_salty_stefan.py` (gate + plots), passing data (31 snapshots), `stefan_salty_{results,timeseries,profiles}.csv`, `fig_*.png`, `run.log`, `analysis_results.txt` |
| `docs/md_files/YANG_ECCO_IMPLEMENTATION_ROADMAP.md` | §A.3b added; A.3/A.2 status annotations; scope + quick-reference updates |
| `docs/md_files/STAGE_A_STEFAN_GATE_REPORT.md` | §6 next-steps updated (intermediate gate done) |

`Boundary.h` remains in the (shared) Stefan-gate configuration; the binary was rebuilt and
verified current against the source tree before the run.

---

## 8. Next steps (per roadmap)

1. **A.2 shakeout** — tiny-grid run with the full Yang flag set (`BOUSSINESQ + EOS_NONLINEAR +
   CONC_VOF_PHASEWEIGHTED + PHASE_CHANGE + ICE_PENALIZATION`, `NConc = 2`), checking NaN-freedom,
   salt confinement to the liquid, budget closure with melting off, and penalization hold. This
   is the only remaining item before Yang.
2. **A.4 Yang production run** — 1440 × 1440 × 4 quasi-2D reference case. The physics gates are
   now closed: any discrepancy against Yang Fig. 2(b)/3(a) is attributable to the *flow-coupled*
   ingredients (EOS, penalization, resolution), not to the melting kernel or the salt coupling.
