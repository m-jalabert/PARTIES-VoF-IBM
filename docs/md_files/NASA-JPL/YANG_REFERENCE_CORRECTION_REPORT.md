# Yang benchmark: reference correction

> **FOLLOW-UP 2026-07-27 (same day): the freshwater "1.78× too fast" identified below is itself a
> RESOLUTION comparison, not a PARTIES defect.** This report's corrected targets stand — figure
> 3(a) is the `ΔSv = 0` series, `t½(Sm=0) = 107.42` — but the conclusion drawn from them in §4/§5
> has been superseded by the convergence ladder run immediately afterwards.
>
> Six rungs (`N` = 216/288/360/512/720/1440) give `t½` = 116.31 / **102.36** / 94.47 / 85.59 /
> 78.89 / 60.29. **Yang's 107.42 falls between the 216 and 288 rungs**, i.e. it equals PARTIES at
> an effective `N = 257` — against Yang's stated **288²** base grid for velocity and temperature.
> Their 5× refinement to 1440² covers only salinity and the phase field, and the `Sm=0` case has
> **no salinity**, so 288² is what sets their heat transport. PARTIES run at that resolution gives
> 102.36 against 107.42, a **4.7 %** difference versus ≈2.8 % digitization uncertainty.
>
> So §5's line *"this is now the largest single discrepancy in the campaign"* should be read as:
> the largest single discrepancy **was** a comparison between PARTIES at 1440² and Yang at an
> effective 288². Full result, mechanism and caveats:
> [YANG_ECCO_IMPLEMENTATION_ROADMAP.md](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) §A.5.2.
>
> **§6's Nusselt anchor does NOT adjudicate — do not cite it as supporting PARTIES.** Inverting
> the lumped closed-box balance for the effective sidewall Nusselt gives:
>
> | | Yang | PARTIES N=288 | N=720 | N=1440 | ladder `t½(∞)` |
> |---|---:|---:|---:|---:|---:|
> | `t½` | 107.42 | 102.36 | 78.89 | 60.29 | ≈63.5 |
> | implied `Nu(θ=1)` | **22.5** | 23.6 | 30.7 | 40.1 | **38.1** |
>
> Churchill–Chu for a vertical wall at `Pr = 10` gives `Nu = 35.6` (`RaT = 1e7`) or `31.4` at the
> effective quadratic-EOS contrast `0.60 RaT`. Converged PARTIES lands on that band and Yang sits
> ~30 % below it — which *looks* like it favours PARTIES, but it does not, because the correlation
> assumes an **unstratified** ambient. Here the box is closed and insulated, and with the density
> maximum at `θ_md = 0.2` the cold meltwater is **denser** than the ambient (`b(0) = −0.04` vs
> `b(1) = −0.64`), so it sinks and builds a stable thermal stratification that suppresses lateral
> transport. The physically correct effective `Nu` should therefore sit *below* Churchill–Chu, and
> a ~30 % suppression is entirely plausible — which is exactly where Yang sits. This anchor is
> consistent with either code being right and settles nothing.

**Controlling finding (2026-07-27): the Stage-A comparison targets were wrong. Figure 3(a) of
Yang et al. (2023) — JFM figure 4(a), the curve the project digitized — is the `ΔSv = 0` series,
not the `ΔSv = 5` series.** The PARTIES production case is `Sm = 5, ΔSv = 5`. It was therefore
compared against a different physical case for the whole campaign, and the freshwater `Sm = 0`
reference was inferred through that same broken chain.

This report supersedes the numerical targets in
[YANG_PRODUCTION_ANALYSIS_REPORT.md](YANG_PRODUCTION_ANALYSIS_REPORT.md),
[YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md), and
`Yang_production/EXPECTED_RESULTS.md`. All PARTIES *measurements* stand; only the reference
values they were judged against change.

It also closes the four questions in [LOHSE_DATA_REQUEST_EMAIL.md](LOHSE_DATA_REQUEST_EMAIL.md)
without author data. That email is no longer blocking.

## 1. Provenance

The arXiv e-print bundle for [arXiv:2302.02357](https://arxiv.org/abs/2302.02357) was retrieved
and archived at `docs/papers/yang2023_arxiv_source/` (`main.tex`, `fig3.pdf`, `fig4.pdf`).
`fig3.pdf` has `sha256 a6bdd059f3fc986ad9b9916e46c163aff6141d6e081731b3b3d9dd7de2c30bd9`.
The figures are **vector** art, so they were rendered locally at 900 dpi and digitized by pixel
tracing — no screen capture, and the only error source is line width. The digitizer is committed
at `Yang_production/digitize_yang_fig3.py` and regenerates every number below.

Built-in calibration check: the `V = 0.5` dashed line of panel (a) maps to the frame mid-height
to `1e-4`.

## 2. What the paper actually says

From `main.tex`, verbatim:

- **`f̄` is defined**: *"we have calculated the average melt rate `f̄ = 1/t_{1/2}`, where `t_{1/2}`
  represents the time needed to melt half of the initial volume"*. **The project's `1/t½`
  convention was correct.**
- **Stratification sense**: the initial condition is
  `S = S_bot + (S_top − S_bot) z/H` in the liquid, and increasing `ΔSv` gives *"stronger stable
  stratification"*. Stable ⇒ salinity increases downward ⇒ `S_bot > S_top`. **PARTIES'
  `cbd2 = 0.5` (top), `cbd5 = 1.5` (bottom) is correct.** (The letter's definition
  `ΔSv = S_top − S_bot` together with positive quoted `ΔSv` values is internally inconsistent;
  the physics and the Huppert–Turner layer-scale agreement both fix the sense unambiguously.)
- **Panel (b) mixes dimensionalities**: *"Circle data points represent 2D simulations, and square
  data points represent 3D simulations."* Panel (a) carries no such statement — which is the
  trap.

## 3. Panel (a) is the `ΔSv = 0` series — three independent proofs

Digitized panel (a) (`Yang_fig3a_reference.csv`):

| `Sm` [g/kg] | `t½` [t_ff] | `V(200)/V₀` | `f̄/f̄₀` |
|---:|---:|---:|---:|
| 0 | **107.42** | 0.2483 | 1.0000 |
| 5 | 118.42 | 0.3142 | 0.9071 |
| 10 | 112.21 | 0.2837 | 0.9573 |
| 15 | 101.84 | 0.2245 | 1.0548 |

Digitized panel (b) 2-D circle series, lightest first (`Yang_fig3b_reference.csv`):

| rank | `Sm=5` | `Sm=7.5` | `Sm=10` | `Sm=15` | identified `ΔSv` |
|---|---:|---:|---:|---:|---|
| 0 (lightest) | **0.910** | 0.939 | **0.955** | **1.048** | 0 |
| 1 | 0.734 | 0.760 | 0.788 | 0.818 | 1.25 |
| 2 | 0.610 | 0.654 | 0.683 | 0.734 | 2.5 |
| 3 | **0.497** | 0.521 | 0.521 | 0.529 | **5** |

1. **The ratios match to 0.007.** Panel (a) gives `f̄/f̄₀ = 0.9071 / 0.9573 / 1.0548` at
   `Sm = 5 / 10 / 15`; panel (b) rank 0 gives `0.910 / 0.955 / 1.048`. No other series is close —
   every other series is at or below 0.82.
2. **Only the `ΔSv = 0` series exceeds 1.0.** Panel (a) has `Sm = 15` melting *faster* than
   `Sm = 0` (`t½` 101.8 vs 107.4). In panel (b) only rank 0 rises above 1.0 at `Sm = 15` (1.048).
3. **The minima land where the paper's own model puts them.** Rank 0's red minimum marker is at
   `Sm ≈ 3.5`, against the weak-stratification prediction `Λ_T = 1` ⇒
   `Sm = Cb ΔT²/(2b₀) = 2.86`. Ranks 3/4/5 have minima at `Sm ≈ 6 / 8 / 9.5`, against the
   strong-stratification prediction `Λ_S = 1` ⇒ `Sm = ΔSv` for `ΔSv = 5 / 7.5 / 10`.

The panel (a) ratios also match the **circles**, not the squares (the 3-D squares on rank 0 read
≈0.845 at `Sm=5` and ≈0.89 at `Sm=10`). **Panel (a) is 2-D**, so PARTIES' `TWOD_CARTESIAN`
comparison is dimensionally correct.

The body text sentence *"In figure 3(a), we plot the normalized volume of ice `V(t)/V₀` … for
different `Sm` with `ΔSv = 5 g/kg`"* is inconsistent with the figure it describes. The figure
wins: it is internally cross-validated against panel (b) three ways.

## 4. Corrected targets and corrected verdict

`f̄₀ = 1/107.42` is the true freshwater normalization, read **directly** off panel (a) rather
than inferred. The `Sm=5, ΔSv=5` target follows from panel (b) rank 3:
`t½ = 107.42 / 0.497 = 216.1`.

| Case | PARTIES | old target | **corrected target** | corrected verdict |
|---|---:|---:|---:|---|
| `Sm=0` (`ΔSv=0`) | `t½ = 60.287` | 59 (inferred) | **107.42** (measured) | **1.78× too fast** |
| `Sm=5, ΔSv=5` | `t½ = 179.848` | 118 | **≈216** | **1.20× too fast** |
| `f̄₅/f̄₀` | 0.3352 | 0.5 | **0.497** | ratio 0.67 of target |
| layer thickness | 0.249–0.296 H | 0.287 H | unchanged | PASS |

The old chain was doubly wrong and the two errors partly cancelled: the `Sm=5, ΔSv=0` curve
(`t½ = 118`) was read as `ΔSv=5`, and the freshwater reference was then *derived* from it as
`0.5 × 118 = 59`. PARTIES' `60.29` matching `59` to 2.2 % was a coincidence of that construction,
not a validation.

**The direction of the error also flips.** PARTIES was believed to melt 52 % *too slowly* in the
salty case. It in fact melts **too fast in both cases**, and much more so in freshwater. The
salty case is now the *closer* of the two.

## 5. What this does to the diagnosis

- **The `f̄₅/f̄₀ = 0.335` vs `0.497` ratio miss is real and unchanged.** Every conclusion built
  on the *ratio* — the τ-null branches, the `Pe_CH`-null branch, the resolved-salt-film audit,
  the Le=100 salt-operator gate — survives untouched.
- **The "freshwater PASS" does not survive.** §A.5.0's exclusion row "Yang geometry + EOS +
  penalization, freshwater — EXCLUDED" is withdrawn. This is now the **largest single
  discrepancy in the campaign** and it involves no salt at all.
- **That is good news for cost.** The freshwater case has `S ≡ 0`, so `Sc = 1000` imposes no
  resolution requirement whatsoever. The `1440²` mesh was demanded by the salt field; the
  freshwater case needs only to resolve `Pe_T = 1e4`. A full freshwater convergence ladder costs
  ≈ 1 300 core-h against ≈ 10 000 for one salty production run.
- **The melting-RB validation is not contradicted, but it is now clearly incomplete.** Favier
  melting-RB (`Pr = 1`, linear EOS, horizontal, `St = 0.1`) passed to 1–3 %. The Yang freshwater
  case differs in three ways at once — `Pr = 10`, the Roquet **quadratic** EOS with the density
  maximum sitting *inside* the temperature range, and lateral (sidewall) rather than
  Rayleigh–Bénard convection. The discrepancy is in that gap.

## 6. Independent order-of-magnitude anchor (does not adjudicate)

A lumped closed-box energy balance — enthalpy `∫θ + F/St` is conserved, so the ambient cools as
ice melts, `θ̄(φ) = (0.5 + 0.4φ)/(1 − 0.1φ)` with `φ = V/V₀` — combined with
`−0.1 St⁻¹ dφ/dt = Nu(θ̄) θ̄ / Pe_T` gives:

| `Nu(θ̄=1)` | 18.6 | 25 | 31.5 | 35.6 |
|---|---:|---:|---:|---:|
| `t½` (laminar `Nu ∝ Ra^{1/4}`) | 130 | 97 | 77 | 68 |

Churchill–Chu for a vertical wall at `Pr = 10` gives `Nu = 35.6` at `RaT = 1e7`, or `31.5` at the
effective quadratic-EOS contrast `0.60 RaT` (`b` spans `−0.04` to `−0.64`, not 0 to 1). That band
(`t½ ≈ 68–77`) sits nearer PARTIES' `60.3` than Yang's `107.4`.

**This does not show Yang is wrong.** The correlation is for an unstratified ambient, and a closed
insulated cavity develops a stable thermal stratification that genuinely suppresses sidewall heat
flux — plausibly by the required factor. It is recorded because it shows the discrepancy lies in a
regime where the modelling assumptions matter, and it is **not** safe to assume PARTIES is the
party in error.

A second measured cross-check points the same way: PARTIES' `max_speed` in the `Sm=0` run is
0.122–0.136 in free-fall units over `t ∈ [4, 60]`, and Yang's figure 4(b) vertical-velocity
profiles peak at ≈0.085. Same order — **the free-fall velocity scale is not off by 1.78**, so the
discrepancy is in heat delivery per unit velocity, not in the nondimensionalization.

## 7. Files

- `docs/papers/yang2023_arxiv_source/` — archived arXiv source (`main.tex`, `fig3.pdf`, `fig4.pdf`)
- `Yang_production/digitize_yang_fig3.py` — reproducible digitizer
- `Yang_production/yang_fig3a_reference.csv` — corrected `V(t)/V₀`, all four `Sm`, `ΔSv = 0`, 2-D
- `Yang_production/yang_fig3b_reference.csv` — `f̄/f̄₀` series
- `Yang_production/yang_fig4_sm5_reference.csv` — **superseded**; it is the `Sm=5, ΔSv=0` curve
  mislabelled as `ΔSv=5`. Retained for provenance only.
