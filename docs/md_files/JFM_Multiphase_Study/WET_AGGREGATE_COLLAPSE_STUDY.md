# Wet Granular Aggregate on a Wall: Simulation Plan

Working title: **Capillary cohesion and collapse of a 2D wet granular aggregate: effect of the contact angle**

Last updated: 2026-09-28

Method: diffuse-interface (Cahn–Hilliard) immersed-boundary method with continuum capillary forcing, as described in
`ARTICLE_JCP_2026/article_V2_07212026/mainV2.tex` (Jalabert & Meiburg, JCP 2026).

Code state assessed: `PARTIES/PARTIES/src` at commit `999c5c5e` (the article code, tag `jcp2026-v2`), plus commit `dc49c334` (2D contact model, which implements C4). This study is developed on branch `VoF-MaximeJalabert`. The phase-change project lives on branch `PhaseChange-ECCO` and is not part of this code state.

Past article inputs used as templates: `/bigscratch/mjalabert314/2D_TESTS/JCP_Article_Sim_2026`.

---

## 1. Summary

We want to simulate 15–20 sand-like disks resting as an aggregate on a wall, in air, with water held in the pore space as liquid bridges. Gravity, liquid–gas surface tension, the solid capillary force (CCF), and grain–grain / grain–wall non-penetration forces are on. There is no cohesive force. The main parameter is the static contact angle `θ`.

**Most of the physics is already in the code.** The planar 2D mode, CH phase field, MCL projection, CCF with torque, gravity, the ACTM + LIN_TAN contact model on the domain walls, and cohesion switched off are all in place. Holding the grains fixed while the bridges form and then releasing them also exists (`STARTUP`).

**Code work needed before production:**

| # | Item | Verdict | Effort |
|---|---|---|---|
| C1 | Air-cavity / liquid-volume drift in closed pockets | **Required.** It affects every run from t = 0, not only after collapse (§6). | ~1 week incl. tests |
| C2 | Initial condition that seeds bridges at every contact of an aggregate | **Required.** No existing `init_type` does this. | 1–2 days |
| C3 | Wall contact angle (bottom wall is currently a 90° wall) | **Decision:** keep a neutral 90° wall (no work) or give the wall the grain `θ` | 3–5 days if needed |
| C4 | Lubrication force in 2D | **Implemented** in `dc49c334` (2D formula × slab, liquid-fraction viscosity); not yet exercised by a test | V2/V4 check only |
| C5 | Runtime switch for mass restoration (hold stage only) | Recommended; moot if C1 option B is done | hours |
| C6 | Diagnostics (per-region volumes, contacts, MCL trace failures) | Recommended | 1 day |

**The key numerical risk is the air-cavity issue.** In 2D, every pore bounded by touching disks is a *closed* air cavity. In hexagonal packing each pore is about 0.03–0.04 D² after bridge seeding. That is exactly the size range where the current scheme removes the higher-pressure phase at a rate of order 25% per 20 gravitational times. The sign of the error flips with `θ`: pores fill artificially with water for wetting grains and bridges drain for non-wetting grains. An uncorrected contact-angle study would therefore be biased in the direction of the trend being measured. §6 gives the mechanism, a quantitative check, and the fix options.

---

## 2. Physical problem

### 2.1 Configuration

- N = 15 or 20 identical rigid disks (sand, `ρ_s/ρ_L = 2.65`) forming a compact aggregate on the bottom wall of a closed box filled with air.
- Water sits in the pore space as pendular bridges at grain–grain and grain–wall contacts.
- Gravity acts downward on grains and fluids.
- Grain–grain and grain–wall contacts use the non-penetration (ACTM normal) and Coulomb friction (LIN_TAN) models. No cohesion model.
- The CCF and IBM pressure load transmit the capillary force to the grains. The liquid–gas surface tension acts in the fluid.
- The main control parameter is the static contact angle `θ` of the grains. Secondary parameters are the Bond number (grain size) and the water content.

### 2.2 What planar 2D implies

These points shape both the physics and the numerics. They are worth agreeing on before running.

1. **Disks are infinite cylinders.** All forces are per unit span. The code integrates them over the storage slab `twod_slab_thickness` (0.01 in the article runs); particle mass and inertia are per slab (`Particle.c`).
2. **A "bridge" at a contact of two touching disks is two separate liquid wedges**, one on each side of the contact point. The contact point disconnects them.
3. **Pores are closed.** Three touching disks in hexagonal packing enclose a pore of area `(√3/4 − π/8) D² = 0.040 D²` (equivalent radius 0.11 D). Four disks in square packing enclose `(1 − π/4) D² = 0.215 D²` (equivalent radius 0.26 D). A pore between two touching bottom grains and the wall is 0.107 D². From the first time step, the aggregate's interior air is a set of small closed cavities. This is why §6 matters for every run.
4. **The 2D meniscus has one curvature, so the capillary pull between touching disks turns repulsive before θ = 90°.** Across the neck plane, the force per wedge is `F' = σ (1 − κ_s y_n)`, where `y_n` is the neck height above the contact point and `κ_s` is the signed meniscus curvature (> 0 convex, liquid pressure above gas). For the equilibrium circular-arc wedge between touching disks, the computed values are:

| filling angle β | θ = 30° | θ = 60° | θ = 90° | θ = 120° | F' = 0 at |
|---|---|---|---|---|---|
| 15° (wedge 0.0013–0.0016 D²) | +21.7 W' | +10.5 W' | −3.6 W' | −16.6 W' | θ = 82.5° |
| 20° (wedge 0.0029–0.0040 D²) | +15.7 W' | +7.0 W' | −3.6 W' | −13.2 W' | θ = 80° |
| 30° (wedge 0.0092–0.0145 D²) | +9.7 W' | +3.6 W' | −3.6 W' | −9.7 W' | θ = 75° |

Here `W'` is the net weight of one grain per unit span for D = 1 mm. Values scale as `1/Bo`, so they are ~9× smaller at D = 3 mm. Numerically the sign change sits at `θ = 90° − β/2`. In 3D, the small-bridge limit between touching spheres is `2πRσ cos θ`, which changes sign at exactly 90°. **Expect the 2D collapse transition between θ ≈ 60° and 90°, not at 90°.** Test V2 checks that the code reproduces this table (§3.4).

### 2.3 Material properties and scaling

| Property | Value |
|---|---|
| Water | ρ = 998 kg/m³, μ = 1.0×10⁻³ Pa·s, σ = 0.0728 N/m |
| Air | ρ = 1.20 kg/m³, μ = 1.81×10⁻⁵ Pa·s |
| Sand (quartz) | ρ_s = 2650 kg/m³ |
| Ratios | `λ_ρ = rho2 = 0.0012`, `λ_μ = mu2 = 0.0181`, `rho_s = 2.655` |

Scaling follows the article's sinking-cylinder case: `L = D`, `U = √(gD)`, `t_g = √(D/g)`. Then `Fr = 1`, `Re` is the Galileo number and `We` is the Bond number `ρ_L g D²/σ`.

| D | Re | We (= Bo) | Oh | t_g | D/ℓ_c | W'/(2σ) | capillary Δt limit at h = 0.01 |
|---|---|---|---|---|---|---|---|
| 0.5 mm | 34.9 | 0.0336 | 0.0053 | 7.1 ms | 0.18 | 0.035 | 5.2×10⁻⁵ |
| **1 mm** | **98.8** | **0.1345** | **0.0037** | **10.1 ms** | **0.37** | **0.14** | **1.0×10⁻⁴** |
| 2 mm | 279.6 | 0.538 | 0.0026 | 14.3 ms | 0.73 | 0.56 | 2.1×10⁻⁴ |
| 3 mm | 513.6 | 1.21 | 0.0021 | 17.5 ms | 1.10 | 1.26 | 3.1×10⁻⁴ |

`W'/(2σ)` is one grain's net weight over the pull of two σ lines. At the proposed baseline D = 1 mm, capillarity dominates by ~7× per contact, so `θ` decides the outcome. D = 3 mm puts gravity and capillarity at the same order and is the natural second Bond number.

---

## 3. Simulation plan

### 3.1 Aggregate geometries

Two candidates. The choice is a decision for §8.

- **A – Pyramid-15.** Hexagonal packing with rows 5-4-3-2-1 (15 disks), bottom row resting on the wall. It is statically stable when dry if base friction suffices. It fails by the base row spreading outward, which the bridges resist. It has 30 grain–grain contacts, 5 grain–wall contacts, 16 closed triangular pores and 4 wall pores.
- **B – Column-20.** Square packing, 4 wide × 5 high (20 disks). Square packing without cohesion is mechanically unstable: it can only stand if the capillary normal forces raise the frictional resistance enough. This is the "collapsing structure" case. It has 31 grain–grain contacts, 4 grain–wall contacts, 12 closed square pores and 3 wall pores.

A perfect lattice under purely vertical gravity can sit in an unstable equilibrium. Give each row of B a small random lateral shift of ≤ 0.01 D, keeping the bottom row fixed. Use **the same seed for every θ**, so that `θ` is the only difference between runs. A ±5% polydispersity is an alternative: `p_mobile.inp` takes a radius per particle.

Generator for `p_mobile.inp` (x, y, z, R per line; z = slab centre as in the article runs):

```python
import numpy as np
R, z = 0.5, 0.005            # D = 1, slab thickness 0.01

def pyramid15(xc=8.0):
    pts = []
    for row, n in enumerate([5, 4, 3, 2, 1]):
        y = R + row * np.sqrt(3.0) * R
        pts += [(xc - (n - 1) * R + 2 * R * i, y) for i in range(n)]
    return pts

def column20(xc=8.0, shift=0.01, seed=1):
    rng = np.random.default_rng(seed)
    s = np.r_[0.0, shift * rng.uniform(-1, 1, 4)]   # bottom row not shifted
    return [(xc + (i - 1.5) * 2 * R + s[j], R + 2 * R * j)
            for j in range(5) for i in range(4)]

def write(pts, fname="p_mobile.inp"):
    with open(fname, "w") as f:
        f.write(f"{len(pts)}\n")
        f.writelines(f"{x:.6f} {y:.6f} {z} {R}\n" for x, y in pts)
```

Grains start exactly in contact with each other and the wall (gap 0). `p_fixed.inp` contains `0`.

### 3.2 Water content

- Seed a fixed liquid area per wedge, `A_w`, at every grain–grain and grain–wall contact. Keep `A_w` **the same for all θ**, so the water content is constant across the sweep.
- Proposed baseline: `A_w = 0.0035 D²`, which corresponds to β ≈ 20° (table in §2.2). Pyramid-15 has about 70 wedges, so ~0.25 D² of water, roughly 20% of the enclosed pore area (pendular regime).
- **Resolution floor:** the contact line on each grain sits at arc length `Rβ` from the contact point, which is 17 cells at β = 20° and D/h = 100. Near the contact point the MCL normals are averaged across the `C_S` ridge between the two grains (normals come from the gradient of the max-combined `C_S`). Keep β ≥ 15–20° so the contact lines stay ≥ ~13 cells away from that ridge.
- Bridges will form only once C2 exists (§5).

### 3.3 Run sequence

1. **Hold, 0 ≤ t ≤ 2** (~5 capillary times `√(ρD³/σ)` = 0.37 t_g). Grains are held at zero velocity (`STARTUP`, `startup_velocity = {0,0,0}`, `startup_time = 2.0`). The seeded liquid relaxes to the θ-equilibrium wedges. Mass restoration is on (static configuration, as in the article's static tests).
2. **Release, 2 < t ≤ 20** (0.18 s physical at D = 1 mm). Grains are free and mass restoration is off, as in the article's moving-body cases, unless C1-B is in place. Collapse of a ~5 D column takes a few free-fall times (`√(2·5) ≈ 3 t_g`). A stable aggregate shows no motion over the remaining time.

### 3.4 Run matrix

**Verification (small domains, short runs; do these first):**

| ID | Purpose | Setup | Pass criterion |
|---|---|---|---|
| V1 | Measure cavity loss (§6.5) | (a) Gas pocket r = 0.05, 0.1, 0.2 D in a liquid box, no grains, no gravity. (b) Triangular pore between 3 fixed touching disks with seeded wedges, θ = 30°, 120°. Run before and after C1. | Pocket volume drift < 1% over t = 20 |
| V2 | 2D bridge force, sign change | Two fixed disks, touching and at gap 0.05 D. One disk on the wall. `A_w` as §3.2, θ = 30, 60, 90, 120° | CCF + IBM load vs `σ(1 − κ_s y_n)` from the measured contour. Sign change near 90° − β/2. Log MCL trace failures in the wedge. |
| V3 | Wall contact angle (only if C3 option ii) | 2D droplet on the flat wall, θ = 30–120° | Angle within ~2% (article: < 1.5% on the cylinder) |
| V4 | Contact mechanics | Dry (no water): one disk dropped on the wall; dry column-20 | Rebound consistent with `e`; overlaps < 1% R; grains at rest do not drift |
| V5 | Hold stage | Chosen aggregate, θ = 60°, t ≤ 2 only | Wedges steady; liquid volume per bridge conserved |

**Production (full domain 16 D × 6 D, D/h = 100):**

| ID | Purpose | Runs |
|---|---|---|
| P-A | Contact-angle sweep (main result) | D = 1 mm; θ = 30°, 60°, 90°, 120°, plus a dry reference. Add θ = 75° once the transition is bracketed (2D prediction: ~80°). |
| P-B | Gravity vs capillarity | D = 3 mm (Bo = 1.21); θ = 30°, 60°, 90° |
| P-C | Water content | θ = 60°, `A_w` × 0.5 and × 2 |
| P-R | Resolution | θ = 60°, D/h = 150 (Cn = 0.75 h), 0 ≤ t ≤ 5 (hold + first 3 t_g after release) |

A D/h = 50 pilot costs ~1/11 of a production run and can screen the geometry choice. The bridges are too coarse at D/h = 50 (meniscus radius ~5 cells), so it is not usable for results.

### 3.5 Measurements

- **Grains** (from `Particle_*.h5`): trajectories, centre-of-mass height, aggregate height and base width (runout) vs time, kinetic energy, and the contact list and coordination number from surface distances.
- **Liquid** (from `Data_*.h5`, connected-component labelling of `C_L ≥ 0.5` and `C_G ≥ 0.5`): number and volume of liquid bodies, bridge rupture and coalescence events, and number and volume of enclosed air pockets. **Per-region volumes are the quality indicator for §6.** The global liquid integral cannot detect the problem: pocket-scale errors are ~0.1% of the global volume in the article's pool cases.
- **Forces per grain:** CCF, hydrodynamic, contact. Check which of these `Particle_*.h5` already stores (C6).
- **Final state:** repose angle (A) or collapse runout (B) vs θ and Bo.

---

## 4. Inputs

### 4.1 Compile-time flags (`src/Include/Boundary.h`)

Relative to `HEAD`:

| Flag | Setting | Note |
|---|---|---|
| `TWOD_CARTESIAN` | define (already) | `AXISYM_RZ` undef |
| No-slip box | already the default branch | Left/right/bottom/top no-slip, `ZPERIODIC` slab. Do **not** reuse the half-domain (free-slip left) builds of the sinking-cylinder runs. |
| `VOF`, `VOF_DIFFUSE`, `SURFACE_TENSION`, `VOF_GRAVITY`, `USE_HYPRE`, `VOF_IBM` | define (already) | |
| `LAG_PARTICLE_RESOLVED`, `SUBSTEP`, `LAG_MARKER_FLAG`, `GRID_UNIFORM` | define (already) | |
| `ACTM`, `LIN_TAN` | define (already) | Normal and tangential contact models; wall collisions active on the x and y walls |
| `COHESION` | undef (already) | No cohesive force, as requested |
| `STARTUP` | **define (change)** | Hold-then-release, §3.3 |
| `LUBRICATION_NORMAL` | define (already) | 2D form since `dc49c334` (§5, C4) |
| `VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE` | 0 for moving grains | Compile-time macro in `VOF_DIFFUSE.c`, not `Boundary.h`; see C5 |

### 4.2 `parties.inp` (baseline: D = 1 mm, θ = 60°, column-20)

```ini
[geometry]
xmin = 0.0
xmax = 16.0
ymin = 0.0
ymax = 6.0
zmin = 0.0
zmax = 0.01

[grid]
NXM = 1600
NYM = 600
NZM = 1
ImportGridFromFile = 0

[twod]
twod_slab_thickness = 0.01
axisym_theta_cells = 1
axisym_theta_span = 6.283185307179586

[simulation]
time_max = 20.0
output_time_interval = 0.02
output_time_interval_2d = 0.02
default_dt = 5.0e-5
max_dt = 5.0e-5
constant_dt = 1
cfl = 0.25
ghost_nodes = 3
resume = 0

[flow]
Re = 98.8
ubulk_target = 0.0
startup_time = 2.0
startup_init = 0
startup_velocity = {0.0, 0.0, 0.0}
vel_init_type = 0

[vof]
We = 0.1345
rho1 = 1.0
rho2 = 0.0012
mu1 = 1.0
mu2 = 0.0181
# new aggregate-bridge initializer, see C2
init_type = <new>
contact_angle_deg = 60.0
weno_order = 5
ch_iter_max = 300
ch_tol = 1.0e-9

[particle]
rho_s = 2.655
grav = {0.0, -1.0, 0.0}
N_forcing_loops = 2
Ndt_coll = 10
e_dry_wall = 0.9
e_dry_particles = 0.9
mu_s = 0.8
mu_k = 0.15
roughness = 3e-3
lub_range = 2.0

[conc]
richardson = {1.0, 0.0}

[lsolve]
P_CG_ETOL = 1.0e-7
P_CG_MAXIT = 300
CG_ETOL = 1.0e-8
CG_MAXIT = 300

[output]
output_vfc = 1
```

Notes:

- `Cn` and `Pe_CH` are left at their defaults (`Input.c`): `Cn = 0.75 h = 0.0075`, `Pe = 0.9/Cn = 120`, as in all article runs. §6 argues for testing a larger `Pe`.
- `Δt = 5×10⁻⁵` is half the capillary limit at D = 1 mm. The adaptive step (`constant_dt = 0`) already includes the capillary constraint `cfl·√((ρ1+ρ2)h³We/4π)` (`Dtime.c`) and is an alternative. At D = 3 mm use `Δt = 1×10⁻⁴`, `Re = 513.6`, `We = 1.21`.
- Since `dc49c334`, `actm_new_pp_collision` uses `e_dry_particles` for grain–grain contacts (it used `e_dry_wall` before). The two are kept equal here anyway.
- `mu_s`/`mu_k` are the code defaults. Whether to use sand-like values (~0.5) is a decision (§8). Hold them fixed across the θ sweep either way.
- Contact stiffness: the ACTM stiffness floor supports one grain weight at 0.1% R overlap (`GRAV_OVERLAP`). Bottom contacts carry several weights plus the capillary pull (up to ~16 W' per wedge at θ = 30°, §2.2). Overlap grows as `(load/W')^{2/3}`, i.e. ~0.5% D for 30 W'. That is acceptable, but log the maximum overlap.

### 4.3 Particle files

`p_mobile.inp` comes from the generator in §3.1. `p_fixed.inp` contains `0`.

### 4.4 Cost

Measured from the article runs on one 40-core node: 0.30 s/step at 307k cells (sinking cylinder), 0.69 s/step at 614k cells (two cylinders), 0.22 s/step at 240k cells (sphere impact). That is **~1.0–1.1 µs per cell per step**.

| Run | Cells | s/step | Steps | Wall time (40 cores) |
|---|---|---|---|---|
| D = 1 mm, 16 D × 6 D, D/h = 100, t = 20 | 0.96 M | ~1.0 | 4.0×10⁵ | ~110 h ≈ 4.6 days |
| D = 3 mm, same grid, Δt = 10⁻⁴ | 0.96 M | ~1.0 | 2.0×10⁵ | ~2.3 days |
| D/h = 150 (P-R), t = 0–5, Δt = 2.8×10⁻⁵ | 2.2 M | ~2.3 | ~1.8×10⁵ | ~4.8 days |
| D/h = 50 pilot, t = 20 | 0.24 M | ~0.25 | ~1.4×10⁵ | ~10 h |

The ~12 production runs (P-A including θ = 75°, P-B, P-C, P-R) total roughly 50 node-days. Runs longer than the 48 h wall limit need `resume`. **Verify that a restart with active contacts reproduces the unrestarted run** (collision history `zeta_t`, contact lists, lagged capillary force); no article run exercised this with grains in contact.

---

## 5. Code changes

### C1 – Air-cavity / liquid-volume drift (required)

See §6 for the mechanism, evidence and options. Recommended path: implement the per-region volume constraint (option B) and use a larger `Pe` (option C) as a sensitivity check.

### C2 – Aggregate bridge initializer (required)

`VoF_Init.c` has initializers for the axisymmetric two-sphere bridge (`VoF_init_axisymmetric_capillary_bridge_two_spheres`), the Nguyen'21 static bridge, and the single- and three-cylinder cases. None of them seeds bridges in an arbitrary packing. Proposed `VoF_init_aggregate_pendular_bridges_2d`:

- Loop over mobile and fixed particles, plus the bottom wall. For every pair with surface distance below a tolerance (e.g. 0.02 D), seed liquid in `{x : max(d_i(x), d_j(x)) < ℓ}`, the points within `ℓ` of both surfaces. For the wall use `d_wall = y − ymin`.
- Solve `ℓ` per contact by bisection so the seeded area per wedge equals a new input `A_w` (`[vof] bridge_area_per_wedge`). The existing `bridge_volume_over_r3` input is specific to the axisymmetric two-sphere case.
- Smooth with the equilibrium tanh profile (width `2√2 Cn`) and clip by `C_S`, as the other initializers do (`vof_init_clip_liquid_by_solid`).
- The seed does not have to be the equilibrium shape: the hold stage (§3.3) relaxes it to the θ-wedge.
- `A_w = 0` gives the dry reference run.

### C3 – Wall contact angle (decision)

Today the bottom wall is a **90° wall**. `VOF_DIFFUSE_set_boundary_values` fills the no-slip-wall ghosts with a zero gradient, and the MCL ghost map (`mcl_global_map_y`) sends ghost rows to row 0. The PLIC wall-wetting path (`VOF_WETTING`) is explicitly blocked with `VOF_DIFFUSE` in `Boundary.h`.

- **(i) Keep the wall neutral (90°).** No work, and it is a clean design for a sweep over grain `θ`: the substrate stays the same. The grain–wall wedges are then asymmetric, with θ on the grain and 90° on the wall.
- **(ii) Wall with its own angle (`contact_angle_wall_deg`, default = grain θ).** Needed if the substrate should be the same material as the grains. Implement a flat-wall wetting ghost fill: the Liu–Ding characteristic construction specialised to a planar wall, or the geometric condition of Ding & Spelt (2007). Apply it consistently in three places: `VOF_DIFFUSE_set_boundary_values` (for `C_L` and `ψ_LG`), the MCL ghost map, and an explicit lagged boundary term in the implicit biharmonic solve. No CCF is needed on the wall because it is fixed. Validate with V3.

Recommendation: start with (i), and revisit after the first P-A results.

### C4 – Lubrication in 2D (implemented in `dc49c334`, to be tested)

**Problem at the article commit.** `lubrication.c` applied the sphere–sphere formula `−6π μ R_eff² g_n / h`, with `μ = 1/Re` (liquid) whatever phase fills the gap. It was not scaled by the slab thickness, while every other force (IBM, CCF, gravity, ACTM stiffness through `M`) is. For two D = 1 disks at gaps of 0.5–2 cells it exceeded the consistent 2D value by 5–10× (slab = 0.01), and the excess scaled with `1/slab`, so results would have depended on the arbitrary storage thickness. No article case had grains in contact, so this was never exercised.

**What `dc49c334` does** (under `TWOD_CARTESIAN`):

- Normal force: the Reynolds-lubrication value for parallel cylinders, `F' = 3√2 π μ g_n (R_eff/h)^{3/2}` per unit span, multiplied by the slab thickness. For a grain–wall contact `R_eff = R`. The gap `h` is floored at `roughness × R`, as in the 3D formula.
- Viscosity: `μ/μ_L = φ + (1 − φ) λ_μ`, where `φ = C_L/(C_L + C_G)` is sampled at the contact point with a 3×3 delta kernel. Solid (`C_S`) tails are excluded from the mix.
- A 2D tangential (Couette) term also exists; it is inactive while `LUBRICATION_TANGENTIAL` stays undefined.

**Still to do:** no test has run this path yet. Check it in V4 (dry disk dropped on the wall: approach and rebound against `e`) and in V2 (two disks in contact inside a bridge). If it misbehaves, the fallback is to undefine `LUBRICATION_NORMAL`: the resolved flow then provides lubrication beyond `lub_range = 2` cells.

### C5 – Runtime switch for mass restoration (recommended)

`VOF_DIFFUSE_RESTORE_MASS_EVERY_STAGE` is a compile-time macro. The two-stage run wants it on during the hold and off after release. Proposed input: `[vof] mass_restore_until = 2.0`. It is not needed if C1-B is implemented, because the per-region constraint then runs throughout.

### C6 – Diagnostics (recommended)

Log at each output: per-region liquid and pocket volumes (reusing the C1-B labelling), contact count and maximum overlap, and MCL trace-failure counters (they exist under `DEBUG_MCL`). Output per-grain CCF, hydrodynamic and contact forces if `Particle_*.h5` does not already store them.

### Not needed

- Gravity on grains and fluid, surface tension, CCF force and torque in planar mode.
- The ACTM + LIN_TAN contact model with domain-wall collisions (x and y walls in the 2D box).
- No cohesion (`COHESION` already undef).
- Per-particle radius in `p_mobile.inp` (for polydispersity).
- Hold-then-release (`STARTUP`, `startup_init = 0`).
- The adaptive time step with the capillary limit.

---

## 6. The air-cavity dissipation problem

### 6.1 Mechanism

At a curved interface, the CH equilibrium has a uniform chemical potential `ψ_0`. The balanced potential-form force converts `ψ_0` into the Laplace jump: `p_L − p_G = (6√2/(Cn We)) ψ_0`. With `p_G − p_L = κ/We` for a gas pocket of curvature `κ`:

```
ψ_0 = −κ Cn / (6√2),   and both bulk values shift by   δ = ψ_0 / F''(bulk) = −κ Cn / (3√2) ≈ −0.24 κ Cn
```

since `F'' = 1/2` at `C_L = 0` and `C_L = 1` for this double well. Inside a gas pocket `C_L` therefore wants to sit at `δ < 0`. Two effects follow.

1. **Dissolution (conservative).** The pocket's gas spreads into the surrounding liquid as a small bulk deficit, so the `C_L = 0.5` contour shrinks even though total gas is conserved. This is the "spontaneous shrinkage" of Yue, Zhou & Feng (2007). Its rate is set by the mobility `1/Pe`.
2. **Clipping sink (non-conservative).** `diffuse_bound_liquid_fraction` clips `C_L` to `[0, 1 − C_S]` after every RK stage. Inside the pocket, the clip resets `δ < 0` to 0 each stage and CH diffusion re-creates it, so gas is destroyed continuously rather than just redistributed.

The same holds with the phases swapped. A convex liquid body (liquid at the higher pressure) sits at `1 + δ`, gets clipped to 1, and loses liquid. **Clipping removes whichever phase is at the higher pressure.** Moving-body runs have mass restoration off (article, *Mass redistribution* subsection), so nothing compensates.

### 6.2 Evidence: 1D radial model with the code's free energy

A 2D-polar Cahn–Hilliard model with the same bulk term, `Cn = 0.0075`, `h = 0.01`, `Pe = 120`, no flow, and a no-flux wall at r = 2 D (script in Appendix A). This is a model of the scheme, not a PARTIES run; V1 measures the real code.

- **Bulk shift check.** The pocket-centre `C_L` settles at −0.0187 and −0.0091, against −0.0207 and −0.0091 predicted by `δ = −κ Cn/(3√2)`.

| Pocket r₀ | Gas left, with clip (as in code), t = 10 / t = 20 | Contour radius at t = 20, clip / no clip |
|---|---|---|
| 0.05 D | pocket gone by t ≈ 10 (also without clip) | 0 / 0 |
| 0.10 D | 0.862 / 0.755 | 0.065 D / 0.086 D |
| 0.20 D | 0.957 / 0.914 | 0.185 D / 0.195 D |

Sensitivity at r₀ = 0.1 D (gas left at t = 20, with clip): baseline 0.755; `Pe × 4` gives 0.928; `h` and `Cn` halved (`Pe = 0.9/Cn`) gives 0.831. A convex *liquid* drop of r₀ = 0.1 D in gas loses exactly the same fraction (0.755), by the `C ↔ 1 − C` symmetry of the double well.

This matches the behaviour already observed: small cavities dragged into the water slowly dissipate. The time unit in the model is the same convective unit as in the runs.

### 6.3 Why it is central for this study

- **Every interior pore is a closed pocket from t = 0 (§2.2).** Hexagonal pores hold ~0.030 D² of air after seeding (equivalent radius ≈ 0.10 D), right in the r₀ = 0.1 D row above.
- **The relevant curvature is that of the corner menisci** (`r_m` = 0.02–0.17 D, §2.2 table), not the pore radius. The shift `δ` is then ~0.01–0.05, so pores can lose faster than the circular estimate suggests.
- **The sign follows θ.** For wetting grains (concave menisci, θ < 90° − β) the pore air is at the higher pressure: air is removed, water is added, and pores fill artificially, so cohesion is **overestimated**. For non-wetting grains (convex wedges) the water is at the higher pressure: bridges drain, so cohesion is **underestimated**. The error amplifies the θ trend being measured.
- **Order of magnitude.** If each of pyramid-15's 16 pores lost ~25% of its air over a run, the added water (~0.12 D²) would be ~50% of the seeded water (~0.25 D²). This estimate is for circular pockets; V1(b) gives the real number.
- **Why the article did not see it.** The article's global conservation figures (≤ 0.03%) cannot detect this: a pocket is ~0.1% of a pool.

### 6.4 Fix options

| Option | What | Removes | Cost / risk |
|---|---|---|---|
| **A. Clip only where needed** | Stop clipping the transported `C_L` to [0, 1] in fluid cells. Keep `C_L ≤ 1 − C_S` in the diffuse solid (needed by MCL/CCF). Clamp `C_L` to [0, 1 − C_S] only when building `ρ`, `μ` and `C_G`. | Clipping sink | Hours. Dissolution remains: in the model the r₀ = 0.1 D contour still shrinks to 0.086 D. Watch every place that assumes `C_L ≥ 0`, especially `C_G = clamp(1 − C_L − C_S)` in the property update: without the clamp, `ρ` could go negative. |
| **B. Per-region volume constraint (recommended)** | Each enclosed air pocket and each liquid body is incompressible and closed, and no evaporation or dissolution is modelled, so **its volume is a physical invariant**. Label connected regions (`C_G ≥ 0.5` and `C_L ≥ 0.5`, excluding `C_S ≥ 0.05`). Track identity across steps by overlap, with the target set at birth; merges sum targets and splits divide by current volume. After each stage, restore each region with the existing `w = C_L(1 − C_L)` redistribution restricted to its own interface band. The open atmosphere is left free. | Clipping sink, dissolution, Ostwald exchange between bridges, MCL-projection loss | ~1 week with tests. Labelling a 2D field of ≤ 1–2 M cells by gathering to one rank is affordable. Main risk: bookkeeping when a grain sweeps the diffuse band next to a pocket (why the article disables global restoration for moving bodies). Here the target is the physical constant rather than a measured integral, but it must be tested (V1 with a moving disk). |
| **C. Lower mobility and/or finer Cn** | Raise `Pe` above `0.9/Cn`, refine `h` and `Cn` | Slows both effects: `Pe × 4` cuts the model loss 24.5% → 7.2%; halving `Cn` cuts it to 16.9% | No code change. Mobility sets how the contact line moves in CH (Yue, Zhou & Feng 2010), and the article's wetting validation used `Pe = 0.9/Cn`, so the droplet-on-cylinder test must be repeated at the new `Pe`. With a capillary-limited `Δt ∝ h^{3/2}`, each halving of `h` costs ~11× in 2D. |

Recommendation: implement **B** (with **A** as a cheap first step that can be tested in a day), and run **C** as a sensitivity check on one production case.

### 6.5 Test V1 (before and after the fix)

- (a) Static gas pocket, r = 0.05, 0.1, 0.2 D, in a liquid box of ~4 D × 4 D, no grains, no gravity, `Pe` = 120 and 480, t = 20. Log the pocket volume (labelled region), `min C_L`, and the per-stage clipped amount.
- (b) Triangular pore between three fixed touching disks with seeded wedges, θ = 30° and 120°. Log pore air volume and wedge volumes: this directly measures the sign-flipping bias of §6.3.
- (c) After C1-B: case (a) with a disk translating past the pocket at the release-stage speeds, to check the moving-solid bookkeeping.

Target: pocket-volume drift < 1% over t = 20.

---

## 7. Other known limitations to keep in view

- **Static contact angle only**, with no hysteresis (article conclusion). Bridge rupture on separation depends on the receding angle in reality.
- **Bridge rupture is diffuse.** The neck pinches when it thins to a few `Cn`, so rupture distance depends on resolution. P-R checks this.
- **Wedge resolution near contacts.** MCL normals are averaged across the `C_S` ridge between touching grains, and near a no-slip wall the MCL samples row 0. Keep β ≥ 15–20° (§3.2) and log trace failures (V2).
- **2D ≠ 3D capillary cohesion** (§2.2). The paper should frame results as planar (liquid prisms between parallel rods), not as spheres.

---

## 8. Decisions needed

1. Geometry: pyramid-15 (A), column-20 (B), or both; perfect lattice with row shifts, or polydisperse.
2. Grain size: D = 1 mm baseline, plus D = 3 mm as the second Bond number?
3. Wall wettability: neutral 90° (no code) or same θ as the grains (C3-ii).
4. Friction and restitution: code defaults (`mu_s = 0.8`, `mu_k = 0.15`, `e = 0.9`) or sand-like values.
5. Water content: single `A_w = 0.0035 D²`, or the P-C sweep.
6. Contact-angle list: 30, 60, 90, 120° plus dry, with 75° added once the 2D transition is bracketed.
7. Cavity fix path: B (recommended) or A + C.

---

## 9. References

- Jalabert, M. & Meiburg, E. (2026). A diffuse-interface immersed-boundary method with continuum capillary forcing for wetting fluid–structure interaction. *J. Comput. Phys.* (manuscript, `mainV2.tex`).
- Liu, H.-R. & Ding, H. (2015); Liu, H.-R. et al. (2017); Washino, K. et al. (2013): as cited in the manuscript (MCL characteristic construction and CCF).
- Yue, P., Zhou, C. & Feng, J. J. (2007). Spontaneous shrinkage of drops and mass conservation in phase-field simulations. *J. Comput. Phys.* 223, 1–9.
- Yue, P., Zhou, C. & Feng, J. J. (2010). Sharp-interface limit of the Cahn–Hilliard model for moving contact lines. *J. Fluid Mech.* 645, 279–294.
- Ding, H. & Spelt, P. D. M. (2007). Wetting condition in diffuse interface simulations of contact line motion. *Phys. Rev. E* 75, 046708.
- Kempe, T. & Fröhlich, J. (2012). Collision modelling for the interface-resolved simulation of spherical particles in viscous fluids. *J. Fluid Mech.* 709, 445–489 (ACTM).
- Biegert, E., Vowinckel, B. & Meiburg, E. (2017). A collision model for grain-resolving simulations of flows over dense, mobile, polydisperse granular sediment beds. *J. Comput. Phys.* 340, 105–127.

---

## Appendix A – 1D radial Cahn–Hilliard model used in §6.2

```python
# 2D-polar CH with the PARTIES free energy:
#   psi = C^3 - 1.5 C^2 + 0.5 C - Cn^2 lap(C),   dC/dt = (1/Pe) lap(psi)
# Gas pocket (C_L = 0) of radius r0 in liquid, no flow, no-flux wall at r = 2.
import numpy as np
h, Cn = 0.01, 0.0075
Pe = 0.9 / Cn
Rdom = 2.0
N = int(Rdom / h)
r = (np.arange(N) + 0.5) * h
rf = np.arange(N + 1) * h
dA = 2 * np.pi * r * h

def lap(q):
    g = np.zeros(N + 1)
    g[1:N] = rf[1:N] * (q[1:] - q[:-1]) / h
    return (g[1:] - g[:-1]) / (r * h)

def run(r0, clip, T=20.0, dt=5e-4):
    C = 0.5 * (1 + np.tanh((r - r0) / (2 * np.sqrt(2) * Cn)))
    gas0 = np.sum((1 - C) * dA)
    for _ in range(int(round(T / dt))):
        psi = C**3 - 1.5 * C**2 + 0.5 * C - Cn**2 * lap(C)
        C = C + dt / Pe * lap(psi)
        if clip:                      # as diffuse_bound_liquid_fraction
            C = np.clip(C, 0.0, 1.0)
    return np.sum((1 - C) * dA) / gas0, C

for r0 in (0.05, 0.10, 0.20):
    for clip in (False, True):
        gas, C = run(r0, clip)
        print(f"r0={r0:.2f} clip={clip}: gas left={gas:.3f}, C_L(centre)={C[0]:+.4f}")
```
