# Yang benchmark — Stage-A validation summary

**Purpose.** The publication-facing summary of the PARTIES validation against Yang et al.,
*Ice melting in salty water*, JFM **969** R2 (2023) / arXiv:2302.02357, prior to the ECCO
sediment-release study. States what was run, what it shows, and what it does not.

**Status 2026-07-30: Stage A CLOSED.** Substantially reproduced, with the limitations of §6 carried
forward — in particular the unexplained −12.8 % salt-stratified melt rate, which closure records
rather than dismisses. Closure unblocks Stage-B implementation and the reduced 3-D shakeout only;
the ECCO production campaign remains gated behind that shakeout. Rationale in the roadmap, §A.5.6.

---

## 1. The reference, and a correction to it

Yang's figures were digitized from the **arXiv vector source** (archived at
`docs/papers/yang2023_arxiv_source/`, `fig3.pdf` sha256 `a6bdd059…`), rendered at 900 dpi and
traced by pixel analysis — not screen capture. Digitizer: `Yang_production/digitize_yang_fig3.py`.
Built-in calibration: the `V = 0.5` dashed line lands on the frame mid-height to 1e-4.

**Figure 3(a) is the `ΔSv = 0` series, not `ΔSv = 5`.** An earlier reading of this campaign took it
as `ΔSv = 5` and derived the freshwater reference from it; both were wrong. Three independent
proofs of the identification are in
[YANG_REFERENCE_CORRECTION_REPORT.md](YANG_REFERENCE_CORRECTION_REPORT.md) §3. The letter's own
body text describes 3(a) as `ΔSv = 5`, which is inconsistent with the figure it labels.

| Yang reference | value |
|---|---:|
| `Sm=0` (`ΔSv=0`) `t½` / `V(200)/V₀` | **107.42** / **0.2483** |
| `Sm=5, ΔSv=0` `t½` / `V(200)/V₀` | 118.42 / 0.3142 |
| `Sm=10, ΔSv=0` | 112.21 / 0.2837 |
| `Sm=15, ΔSv=0` | 101.84 / 0.2245 |
| **`Sm=5, ΔSv=5`** (the production case) | **`f̄/f̄₀ = 0.497` ⇒ `t½ = 216 ± 11`** |

`f̄ = 1/t½` is confirmed verbatim from the paper source. Note the asymmetry in reference quality:
the freshwater targets are **directly digitized curves** (±3, 2.8 %); the `ΔSv=5` target is
**derived** from one marker and one curve, carrying **±11 (4.9 %)**.

## 2. The decisive methodological fact

**Yang uses two grids; PARTIES uses one.** Their `u,T` sit on **288²**; only `S` and `φ` are
refined 5× to 1440². PARTIES is single-grid, so no PARTIES run can match both at once — it is
either over-resolved in `u,T` (at 1440²) or under-resolved in `S,φ` (at 288²).

This matters because the two cases behave oppositely:

| | freshwater | salty (`ΔSv=5`) |
|---|---|---|
| `t½` across `N` = 288 → 1440 | 102.36 → 60.29 (**70 %**) | 188.39 → 179.85 (**5.4 %**) |
| verdict | strongly resolution-sensitive | **converged** |

Physical reading: `ΔSv = 5` stratification suppresses convection (the paper's own mechanism for the
melt-rate minimum), so the salty flow is laminar enough to be resolved at 288²; the freshwater case
has vigorous convection whose heat delivery is resolution-limited.

## 3. Runs published

All on binary `6a8fa6eb` unless noted; provenance chain in §5.

**(a) Freshwater convergence ladder** — the central figure.
`ECCO_TESTS/StageA_Yang/Yang_fresh_convergence/`, analyser `analyze_ladder.py`.

| `N` | 216 | 288 | 360 | 512 | 720 | 1440 |
|---|---:|---:|---:|---:|---:|---:|
| `t½` | 116.31 | **102.36** | 94.47 | 85.59 | 78.89 | 60.29 |

Monotone and predictive: it forecast `N720` at 79.5 (got 78.89, 0.8 %) and `N288` at 101–102 (got
102.36). **Yang's 107.42 falls between the 216 and 288 rungs — effective `N = 257`, against their
stated 288² base grid (10.8 %).**

**(b) Salty resolution ladder** — establishes that (a)'s argument does *not* transfer.
`ECCO_TESTS/StageA_Yang/Yang_salty_resolution/`, analyser `analyze_salty.py`.
`t½` = **188.39 / 178.28** at `N` = 288 / 432 (both `6a8fa6eb`), i.e. **5.4 % across a 1.5×
refinement**, against 12 % for freshwater over the same interval. The 1440² value 179.85 is
consistent but its binary provenance is unverifiable (§5), so it corroborates rather than counts.

**(c) Matched-resolution comparison at 288²** — the validation claim.

| | PARTIES 288² | Yang | difference |
|---|---:|---:|---:|
| fresh `t½` | 102.36 | 107.42 | **−4.7 %** |
| fresh `V(200)/V₀` | **0.2447 ± 0.019** | 0.2483 | **−1.4 %** |
| salty `t½` | 188.39 | 216 ± 11 | **−12.8 % (2.6σ)** |
| `f̄₅/f̄₀` | 0.543 | 0.497 | **+9.3 %** |
| layer spacing | 0.215–0.380 H | 0.287 H | brackets (see §6) |

**(d) Converged PARTIES answer at 1440²** — what the method gives fully resolved.
Freshwater `t½ = 60.29`, `V(200)/V₀ = 0.1979`; salty `t½ = 179.85`.

## 4. Supporting verification (previously passed, unchanged)

1-D Stefan (`λ` to 1.7 %); 1-D binary/salty Stefan; nonlinear-EOS reconstruction; Brinkman
penalization sweep; Le=100 salt-operator gate against an exact two-phase similarity solution;
Favier melting-RB at `St=0.1` (plateau ratio 0.98–1.03, budget `1.0006 ± 0.013`); flat-wall 2-D RB
(`Nu = 8.27 ± 0.26` at `Ra=1e6`).

## 5. Binary provenance

Long runs spanned three builds, so equivalence was measured, not assumed:

| link | evidence |
|---|---|
| `aac878bf` (1440² freshwater production) ≡ `6a8fa6eb` | `Yang_fresh_provenance/`: **1.7e-12** at `t=2`, 5.8e-5 at `t=6`, vs a 7.2e-3 discriminating size — **125× margin**; also clears rank-count independence (256 → 128) |
| `675f4f56` (07-18) vs `6a8fa6eb`, **salty** | **TEST INVALID — not a lineage result.** `675f4f56` is `parties.melting_rb_backup`, compiled with **periodic x-boundaries** for the Favier melting-RB campaign; the Yang case needs no-slip walls. The 9.2 % `V` difference it produced measures that BC mismatch, not the binaries. See note below |
| `6a8fa6eb` ≡ `e3ec86e2` (flag off) | gate bisect: identical to every printed digit |

**The salty 1440² production run's binary lineage is UNVERIFIABLE with what survives.** The binary
that produced it no longer exists, its trajectory is on an unreachable filesystem, the source
predates any commit (the Yang work is entirely uncommitted), and the one surviving older binary is
compiled for a different problem. No test can close this link.

**This does not affect the published claims, because none of them load-bear on that run.** The
salty comparison point is `N288` and the salty convergence claim rests on `N288` + `N432`, all on
`6a8fa6eb`. The 1440² salty `t½ = 179.85` is reported as **corroborating, provenance unverified**,
not as evidence. Consequence to state honestly: the salty convergence claim is **two verified
points** (188.39 → 178.28, 5.4 % across a 1.5× refinement) plus one unverified third, not three
points across a 5× range.

## 5b. Carry-over to the Stage-B build

The Stage-A results above were produced with `VOF_IBM` and `LAG_PARTICLE_RESOLVED` **off**. Every
Stage-B run has them **on**, so their non-invasiveness was measured rather than assumed
(roadmap §B.0, `ECCO_TESTS/StageB1_development/StageB_noninvasive/`):

- **Le=100 salt gate** (quiescent): `λ = 0.21844955`, profile `L∞ = 1.441e-2`, salt drift
  `2.764e-3` — identical to the Stage-A reference to every printed digit.
- **Freshwater 288² to `t=10`** (flow-coupled, `max|u| = 0.2009`): worst relative ice-volume
  difference **8.9e-9**, enthalpy drift identical to four digits. Six orders of magnitude tighter
  than the accepted Stage-A binary-lineage test at the same time.

The build is not bit-identical to Stage A and cannot be (§5), but it is indistinguishable from a
rebuild of the same physics by the same measure this campaign used to certify its own lineage.

## 6. Limitations — to be stated in the paper, not resolved by more compute

1. **The salty absolute misses the project's own 10 % acceptance**: −12.8 % at matched resolution,
   2.6σ against a target that is itself derived and carries ±4.9 %. The single-grid argument does
   **not** cover this — the salty case is resolution-converged, so the grid split cannot be blamed.
   This is a genuine unexplained residual.
2. **The 288² agreement may be mechanistic or coincidental.** PARTIES at 288² carries a diffuse
   band **5.3× wider** than Yang's (§G.8 of the roadmap: 5–95 % band 2.17e-2 vs 4.09e-3), so both
   codes under-deliver interfacial heat at that effective resolution, plausibly for different
   reasons. The experiment designed to separate this (roadmap §A.5.5 item 3) failed on an enthalpy
   leak and has not been redone.
3. **Morphology at 288² is a consistency check, not a like-for-like comparison** — only 3–4 layers
   are resolvable and the estimator used differs from the one behind the 1440² figures.
3b. **`V(200)/V₀` carries a ±0.019 systematic from late-time enthalpy drift.** The `t½` measurement
   is clean (drift 2.0e-5 at `t=102`), but drift grows after the ice melts through to the wall:
   2.3e-4 at `t=140`, 8.8e-4 at `t=170`, **1.68e-3 at `t=200`**. Since `E = ∫θ + ∫F/St`, that last
   figure could account for `±0.019` in `V/V₀` — i.e. **±7.7 % of the reported 0.2447**. The −1.4 %
   agreement with Yang is therefore *within* the run's own systematic, and should be quoted as
   consistent rather than as a 1.4 % measurement. Closing this needs the late wall-melt-through
   audit (roadmap §A.5.1/§7.4), which remains open.
4. **Band-localized model corrections cannot fix (1).** The Brinkman drag `ρ(2/τ)φ_s` overlaps the
   diffuse band (1e4 at `F=0.5`, 2e3 at `F=0.9`), so buoyancy changes there are damped. Verified by
   the variant-1 test: a change acting only in the band moved `t½` by **−0.09** against the +37.8
   needed. Yang's `η = dt ≈ 1e-4` is the same strength, so this constrains **both** codes.

## 7. What is *not* claimed

- That the 1440² salty production run reproduces Yang. It compares PARTIES-converged against
  Yang-at-288² and is reported as the converged answer, not as the validation.
- That Yang's published rates are wrong. They are consistent with their stated base-grid resolution.
- Any converged `t½` extrapolation. The ladder is not converged anywhere in the sampled range, and
  the earlier `t½(∞) ≈ 63–64` figure came from a fit extrapolated a factor of two beyond its data.
