# Yang finite-interface salt-operator gate

> **§7 SUPERSEDED 2026-07-27 — the recommendation, not the evidence.** The project goal is to
> reproduce Yang **with PARTIES**; an AFiD-MuRPhFi run is not a deliverable. Section 7's
> "Recommended route: reproduce Yang with AFiD-MuRPhFi" is **withdrawn**, and with it the framing
> of a PARTIES result as conditional on a prior Allen–Cahn port. Sections 1–6 stand unchanged and
> remain the controlling evidence: the salt-only graft is rejected, and note that the PARTIES
> Cahn–Hilliard formulation **beat** the mapped Yang operator against the exact similarity
> solution. The replacement plan is
> [YANG_ECCO_IMPLEMENTATION_ROADMAP.md](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) §A.5, which tests
> grid/interface-width convergence and the `f̄(Sm)` curve shape before any interface-model rewrite.

**Controlling decision (2026-07-21): STOP. Do not run the refined 1-D gate,
the proposed `Sc=500` branch, another PARTIES Yang production case, or a
parameter sweep.** The allocation-free matched-state audit already rejected an
under-resolved salt film, and the two cheap 1-D tests in this report reject a
salt-equation-only correction. An accurate Yang reproduction now requires the
**coupled Allen--Cahn/Hester phase and salinity method**, not another adjustment
to `Pe_S`, Darcy damping, Cahn--Hilliard mobility, or the salt source alone.

This report is the controlling continuation of
[YANG_PRODUCTION_ANALYSIS_REPORT.md](YANG_PRODUCTION_ANALYSIS_REPORT.md) and
[YANG_SM0_ANALYSIS_REPORT.md](YANG_SM0_ANALYSIS_REPORT.md). It records the
smallest additional calculation justified by the transferred
`Results_from_POD` products and closes the salt-only branch before it consumes
production allocation.

## 1. Evidence entering this gate

The completed `Sm=5` production run melted too slowly:

- PARTIES `t_1/2 = 179.8484`, compared with the digitized Yang target `~118`;
- the matched fresh run passed, with `t_1/2,0 = 60.2866` versus the inferred
  Yang value `~59`;
- consequently PARTIES gives the salty/fresh melt-rate ratio `0.3352`, versus
  Yang's plotted value of approximately `0.5`.

The data-local matched-volume audit then established that the discrepancy is
not an unresolved numerical salt film:

- snapshot source and finite-difference melt rate agree within 0.53 %;
- `Sm5/Sm0` melt-rate ratios are 0.25--0.37, heat-flux ratios are 0.40--0.47,
  and interface-superheat ratios are 0.25--0.36;
- all salinity profiles are usable, none has `delta_S,90 < 4 dx`, and median
  recovery lengths are 27.7--98.4 cells;
- salt overlap/leakage indices are only `1.25e-4`--`1.92e-4`, with no starved
  or refreezing rows under the predeclared definitions.

That audit rejected the conditional `Sc=500` calculation. It instead motivated
one cheap check of the remaining equation-level difference: Yang's
finite-interface salinity operator versus the legacy PARTIES masked-diffusion
operator.

## 2. Source audit and model distinction

The reference implementation audited here is the official
[AFiD-MuRPhFi repository](https://github.com/chowland/AFiD-MuRPhFi), commit
`f5389960543d7e379a9d7562ccee4bfa23fd0a8c`. With Yang's solid phase `phi` and
PARTIES liquid fraction `F = 1-phi`, its salinity equation is

```text
dS/dt + u.grad(S)
  = (1/Pe_S) lap(S)
    + [(1/Pe_S) grad(F).grad(S) - S dF/dt] / (F + delta),
delta = 1e-6.
```

Equivalently, the two diffusion terms can be written

```text
(1/Pe_S)/(F+delta) div[(F+delta) grad(S)].
```

AFiD evaluates explicit salinity terms, advances the Allen--Cahn phase field,
and then obtains the latent-salt term from the **complete implicit phase
increment** before its implicit salinity solve. Its published multicomponent
validation starts from nonzero-time analytic temperature and salinity profiles;
the numerical salinity is continuous across the diffuse boundary even though
the physical solid concentration is interpreted as zero. See the project's
[phase-field validation](https://chowland.github.io/AFiD-MuRPhFi/phasefield_validation/).

Legacy PARTIES instead evolves a conservative Cahn--Hilliard liquid fraction
and diffuses salt with a masked flux proportional to `F grad(S)`. Its physical
melt source is distinct from phase-band relaxation and boundedness projection.
Therefore replacing the salt equation alone cannot reproduce AFiD's coupled
discrete phase/salt update exactly.

## 3. Experimental implementation

An optional, default-off runtime path was added behind:

```ini
[phase_change]
yang_salt_transport = 0
yang_salt_delta = 1.0e-6
```

When enabled for concentration field 1, it:

1. evaluates the combined diffusion in the stable flux form above;
2. includes that operator consistently in both Crank--Nicolson halves;
3. uses the `(F+delta)`-weighted inner product that makes the row-scaled
   diffusion matrix self-adjoint for the existing CG solver; and
4. adds `-S*melt_src/(F+delta)` with the existing low-storage RK history.

Only the physical melt source is used for salt rejection. Treating the full
Cahn--Hilliard phase increment as latent melt caused an immediate nonphysical
interface mode because CH relaxation and admissibility correction do not
consume latent heat. Splitting `grad(F).grad(S)/(F+delta)` explicitly was also
unstable on the clipped CH band. Those two diagnostic attempts were stopped
near `t = 0.006`, and their small outputs were retained as provenance.

The stable flux-form code built cleanly and leaves the production default off.
It should be described as a **Yang salt operator mapped onto PARTIES**, not as
the complete Yang method.

## 4. Gate design

The gate is under
`PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_salt_Le100_gate/`. It uses a quiescent
`768 x 4 x 1` strip with

```text
Pe_T = 100, Pe_S = 10000, Le = 100, St = 0.5,
liquidus slope = 0.5, delta = 1e-6, t_end = 0.5.
```

The exact two-phase similarity problem gives

```text
lambda          = 0.21713826597063368
S_interface     = 0.0011640201641615062
theta_interface = -0.0005820100820807531.
```

Each comparison changes only `yang_salt_transport`:

- **production-like start:** a diffuse phase profile with salt initialized
  consistently with the initial PARTIES liquid mask;
- **smooth similarity start:** both branches restart at `t=0.25` from the same
  analytic temperature and continuous salinity profiles and a diffuse front at
  the exact similarity position, then advance to `t=0.5`.

The second test removes the initial salt jump as a confounder. The predeclared
direction rule allowed a refined test only if the mapped branch remained
finite, conserved physical salt within 1 %, and improved either the final
liquid-salinity `L_inf` error or the fitted front constant.

## 5. Results

### 5.1 Production-like initialization

Slurm job `19427758` completed in 2 min 15 s on four ranks.

| Metric | legacy | mapped Yang salt | finding |
|---|---:|---:|---|
| fitted `lambda` | 0.21844955 | 0.21908046 | both close |
| relative `lambda` error | 0.6039 % | 0.8945 % | mapped branch 1.48x worse |
| final liquid-S `L_inf` | 1.441e-2 | 8.773e-2 | mapped branch 6.09x worse |
| maximum fitted-window S error | 2.781e-2 | 1.215e-1 | mapped branch worse |
| physical-salt drift | 2.764e-3 | 2.599e-4 | mapped branch 10.6x better |
| enthalpy drift | 4.98e-10 | 4.78e-10 | both negligible |
| maximum speed | 0 | 0 | quiescence exact |

The conservative quantity improves, but neither accuracy target improves.

### 5.2 Smooth nonzero-time similarity initialization

Slurm job `19427886` completed in 2 min 20 s on four ranks.

| Metric | legacy | mapped Yang salt | finding |
|---|---:|---:|---|
| fitted `lambda` | 0.21838506 | 0.21849378 | both close |
| relative `lambda` error | 0.5742 % | 0.6243 % | mapped branch 8.7 % worse |
| final liquid-S `L_inf` | 4.371e-4 | 8.719e-4 | mapped branch 1.99x worse |
| physical-salt drift | 9.458e-6 | 1.356e-5 | mapped branch 1.43x worse |
| plain-S drift | 2.44e-15 | 1.898e-5 | legacy alone conserves plain S |
| enthalpy drift | 1.90e-10 | 1.89e-10 | both negligible |
| maximum speed | 0 | 0 | quiescence exact |

Both branches are finite; the mapped branch salinity remains in
`[4.18e-5, 1]`. Removing the startup discontinuity therefore does not rescue
the salt-only mapping. The difference is small in the front but systematic and
is paired with a factor-of-two profile penalty.

### 5.3 Allocation used

All jobs, including the two early-stop diagnostics and one module-load failure,
consumed about **0.48 core-hour total**. The two completed scientific gates
consumed `9:00 + 9:20 = 18:20` CPU time, or 0.306 core-hour. The failed explicit
forms were cancelled at their first instability, and the pending duplicate
used no allocation.

## 6. Interpretation

The result separates three claims:

1. **The production anomaly is real.** The freshwater control, salty/fresh
   ratio, and matched transport audit all agree that PARTIES supplies too little
   convective heat to the salty interface compared with the digitized Yang
   result.
2. **It is not a salt-resolution failure.** The production salinity boundary
   layer has tens of cells, leakage is tiny, and `Sc=500` is rejected.
3. **A salt-only equation graft is not the correction.** It is numerically
   stable in flux form, but loses to legacy PARTIES from both a production-like
   and an analytic smooth start. The residual method difference is the coupled
   phase evolution and its substep-level salt/latent-heat update, principally
   Allen--Cahn/Hester in AFiD versus conservative Cahn--Hilliard in PARTIES.

This does not prove that the complete Yang/AFiD formulation is wrong or that
the published melt curve is exact. It proves that enabling the optional
salt-only PARTIES path cannot be presented as an accurate reproduction and
does not justify a production run.

## 7. Allocation-minimizing path to an accurate Yang simulation — WITHDRAWN 2026-07-27

> **This entire section is superseded** by
> [YANG_ECCO_IMPLEMENTATION_ROADMAP.md](YANG_ECCO_IMPLEMENTATION_ROADMAP.md) §A.5. It is retained
> verbatim below only as provenance for the 2026-07-21 decision. Do **not** act on it: the
> benchmark must be reproduced with PARTIES, and §A.5 orders the remaining hypotheses by evidence
> per core-hour, with the Allen–Cahn port ranked last rather than first. The `t=60` direction test
> and its 25 %-of-gap threshold defined below are carried over into §A.5 and remain in force.

### ~~Recommended route: reproduce Yang with AFiD-MuRPhFi~~ (withdrawn — off-goal)

If the immediate goal is an accurate **Yang reproduction**, use the official
AFiD-MuRPhFi method rather than continuing to hybridize PARTIES. This preserves
the paper's Allen--Cahn phase equation, full phase increment in latent salt and
heat, substep order, phase-grid/salinity-grid refinement, and direct solid
forcing. Pin the audited commit above and first reproduce its bundled 1-D
multicomponent validation. This is the lowest-risk route because it tests the
reference implementation rather than a new cross-code discretization.

The compute gate is:

1. perform build/input/schema checks without allocation;
2. run **one coarse bundled 1-D validation**;
3. only if its liquid temperature, liquid salinity, and `sqrt(t)` front errors
   match the repository validation, run **one refined 1-D confirmation**;
4. only after both pass, launch one 2-D Yang trajectory initially ending at
   `t=60`; if the predeclared `t=60` direction test passes, extend the same
   checkpointed trajectory to `t=200` rather than starting another case.

For the `t=60` direction test, compare against the existing PARTIES value
`V/V0=0.77427` and digitized Yang value `0.69664`. Require the new method to
close at least 25 % of that gap (`V/V0 <= 0.75486`) while maintaining bounded
phase/salinity, solver convergence, salt and enthalpy closure, and the expected
layer morphology. Stop at `t=60` if it does not. Passing this direction gate is
permission to extend, not final benchmark acceptance; final acceptance still
requires the `t_1/2`, `V(200)/V0`, and layer-spacing criteria.

### PARTIES-only route: port the coupled method before any new physics run

If the scientific requirement is specifically to reproduce Yang **inside
PARTIES**, first implement the Allen--Cahn/Hester phase equation and couple its
complete per-stage phase increment to latent heat and salt in the AFiD order.
The implementation must include the same liquidus forcing, interface salt-flux
term, `delta`, phase thickness/mobility mapping, and solid-velocity treatment.
Do not combine the present experimental salt switch with the production CH
model and call it Yang.

Before allocating a job, verify the discrete operators on manufactured arrays,
phase/salt sign conventions, single-rank versus four-rank identity, restart
history, and conservation identities. Then apply the same one-coarse/one-refined
1-D gate and the single checkpoint-extendable `t=60` 2-D direction run above.

This route costs more development time and carries more verification risk than
running AFiD, but it is the only defensible way to claim a Yang-method result
from PARTIES.

## 8. Files and provenance

Primary gate artifacts:

- `PARTIES/testcases/ECCO_TESTS/StageA_Yang/Yang_salt_Le100_gate/gate_comparison.csv`
  (`SHA-256 9d008445ec108b5dc40dd7420f70349fc83d04b372da4371aee6284367cc882a`);
- `similarity_gate_comparison.csv`
  (`b2286056048faf40450b67fa31ef74549cd531d6025949a2c5899b8b0b3245a7`);
- `analyze_gate.py`
  (`5079cd51625eb2b67314f35840782c5a667e9e53d0478fbc83001ccde2479f3f`);
- `make_similarity_seed.py`
  (`122eb62cbd23f9f4db60126446dbd4bdc22966609408e79509ba44520c317d87`).

Scientific references:

- [Yang et al., JFM 969 R2 (2023)](https://doi.org/10.1017/jfm.2023.582);
- [arXiv manuscript](https://arxiv.org/abs/2302.02357);
- [AFiD-MuRPhFi source](https://github.com/chowland/AFiD-MuRPhFi).
