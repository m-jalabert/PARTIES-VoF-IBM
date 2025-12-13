# VOF-IBM Architecture Description

This document outlines the code architecture and physical methods used when the simulation is configured for **Volume of Fluid (VOF)** coupled with **Immersed Boundary Method (IBM)**, including **Surface Tension** and **Moving Contact Lines**.

**Configuration Flags:**
- `VOF_PLIC`: Geometric VOF with Piecewise Linear Interface Calculation.
- `SURFACE_TENSION`: Surface tension enabled (CSF method).
- `VOF_GRAVITY`: Gravity enabled for two-phase flows.
- `VOF_IBM`: VOF-IBM coupling for wetting on immersed solid boundaries.
- `LAG_PARTICLE_RESOLVED`: Resolved Lagrangian particles (acting as IBM solids).

---

## 1. Initialization Phase (`Cart3d.c`)

The initialization process sets up the Eulerian grid, Lagrangian particles, and VOF fields.

### 1.1 Memory Allocation

- **VOF Fields**: `VoF_create()` allocates:
  - Volume fraction $F$ and smoothed version $\tilde{F}$
  - Interface normals $(\mathbf{n}_x, \mathbf{n}_y, \mathbf{n}_z)$ and smoothed normals
  - Plane intercept $\alpha$, curvature $\kappa$
  - Fluxes, convective terms, density $\tilde{\rho}$, viscosity $\tilde{\mu}$
  
- **IBM Fields** (when `VOF_IBM` enabled):
  - Solid volume fraction $vfc$ and smoothed version $vfc_{smooth}$
  - Solid normals $(\mathbf{n}_{s,x}, \mathbf{n}_{s,y}, \mathbf{n}_{s,z})$
  - CCF tangent vectors $\mathbf{t}_{int}$, $\mathbf{t}_{cl}$
  - Contact line force density $\mathbf{f}_{CCF}$

### 1.2 Particle Initialization

`Particle_initialize()` sets up resolved particles defining the solid domain:
- Reads particle positions, radii, and properties from input files
- Generates Lagrangian surface markers on each particle
- Each marker has:
  - Position $(X_L, Y_L, Z_L)$
  - Associated volume $V_L = 4\pi R^2 / N_L$ (surface area per marker)
  - Flag for active/inactive status

---

## 2. Immersed Boundary Method (IBM)

The IBM enforces no-slip boundary conditions on moving solid boundaries using a **direct forcing approach** with Lagrangian surface markers.

### 2.1 Delta Function Interpolation

The Eulerian and Lagrangian representations are coupled via a regularized delta function. The code uses the **Roma kernel** (Roma et al., 1999):

$$\delta_h(r) = \begin{cases} 
\frac{1}{3}\left(1 + \sqrt{1 - 3r^2}\right) & |r| \leq 0.5 \\
\frac{1}{6}\left(5 - 3|r| - \sqrt{1 - 3(1-|r|)^2}\right) & 0.5 < |r| \leq 1.5 \\
0 & |r| > 1.5
\end{cases}$$

The 3D kernel is constructed as a tensor product:
$$\delta_h(\mathbf{r}) = \delta_h(r_x/h) \cdot \delta_h(r_y/h) \cdot \delta_h(r_z/h)$$

This kernel has a support width of $3h$ and satisfies the moment conditions for 2nd-order accuracy.

### 2.2 Interpolation: Eulerian → Lagrangian (`Interpolate_Eul_to_Lag`)

Interpolates an Eulerian field $a(\mathbf{x})$ to Lagrangian marker positions:
$$A_m = \sum_{i,j,k} a_{i,j,k} \cdot \delta_h(\mathbf{X}_m - \mathbf{x}_{i,j,k})$$

For each velocity component, the interpolation uses the appropriate staggered grid locations (xu for u, yv for v, zw for w).

### 2.3 Spreading: Lagrangian → Eulerian (`Interpolate_Lag_to_Eul`)

Spreads Lagrangian quantities back to the Eulerian grid:
$$a_{i,j,k} = \sum_{m} A_m \cdot \delta_h(\mathbf{X}_m - \mathbf{x}_{i,j,k}) \cdot \frac{V_L}{h^3}$$

where $V_L/h^3$ is the dimensionless marker volume.

### 2.4 IBM Forcing (`Lagrangian_force`, `Lagrangian_force_individual`)

The IBM forcing is applied in a predictor-corrector framework within each RK3 substep:

**Step 1: Predictor** (before implicit viscous solve)
- Called with `Lagrangian_force(-1, ...)` to apply initial forcing
- Solves explicit velocity update first
- Computes IBM force and adds to RHS

**Step 2: Corrector** (after implicit viscous solve)
- Called with `Lagrangian_force(N_forcing_loops, ...)` for iterative correction
- Iterates to enforce no-slip more accurately

**For each particle and each velocity component:**

1. **Interpolate velocity to markers:**
   $$U_{L,m} = \sum_{i,j,k} u_{i,j,k} \cdot \delta_h(\mathbf{X}_m - \mathbf{x}_{i,j,k})$$

2. **Compute desired marker velocity** (rigid body motion):
   $$U_d = \mathbf{U}_p + \boldsymbol{\Omega}_p \times \mathbf{r}_m$$
   where $\mathbf{r}_m = \mathbf{X}_m - \mathbf{X}_p$ is the position relative to particle center.

3. **Compute Lagrangian acceleration:**
   $$F_L = \frac{U_d - U_{L,m}}{2\beta\Delta t}$$
   
   where $\beta$ is the RK3 coefficient for the current substep.

4. **Accumulate reaction force on particle** (for mobile particles):
   $$\mathbf{F}_p -= \rho_L \cdot M_L \cdot F_L$$
   $$\mathbf{T}_p += \rho_L \cdot M_L \cdot (\mathbf{r}_m \times F_L)$$
   
   where $\rho_L$ is the local fluid density (from `VOF_IBM` coupling) and $M_L$ is the marker mass.

5. **Spread force to Eulerian grid:**
   $$f_{i,j,k} = \sum_m F_L \cdot \delta_h(\mathbf{X}_m - \mathbf{x}_{i,j,k}) \cdot \frac{V_L}{h^3}$$

6. **Add to momentum RHS with density weighting** (when `VOF_IBM`):
   $$\text{RHS}_u += 2 \rho_{face} \cdot f_{i,j,k}$$
   
   where $\rho_{face}$ is the face-centered mixture density.

### 2.5 Solid Volume Fraction (`Interpolate_add_to_volume_fraction`)

Computes the solid volume fraction $\Phi_s$ at each grid location:
1. For each particle, find bounding box
2. Compute level set function $\phi = |\mathbf{x} - \mathbf{X}_p|/R - 1$
3. Use sub-cell integration to determine volume fraction:
   $$vf_{i,j,k} = \frac{\sum_{\phi < 0} |\phi|}{\sum |\phi|}$$

This is used for density/viscosity blending in `VOF_IBM` mode.

### 2.6 Momentum Integration (`Interpolate_integrate_momentum`)

Computes the integrated fluid momentum over the particle volume:
$$\int_{\Omega_p} \tilde{\rho} \mathbf{u} \, dV = \sum_{i,j,k} \rho_{face} \cdot u_{i,j,k} \cdot \delta_h(\mathbf{X}_m - \mathbf{x}_{i,j,k}) \cdot h^3$$

The face-centered density is computed as:
$$\rho_{face,u} = \frac{1}{2}(\rho_{i,j,k} + \rho_{i-1,j,k})$$

---

## 3. Time Integration Loop (`Temporal_int.c`)

The simulation uses a **3-stage Runge-Kutta (RK3)** scheme with coefficients:

$$\gamma = \left(\frac{8}{15}, \frac{5}{12}, \frac{3}{4}\right), \quad \zeta = \left(0, -\frac{17}{60}, -\frac{5}{12}\right), \quad \beta = \left(\frac{4}{15}, \frac{1}{15}, \frac{1}{6}\right)$$

Each RK substep calls `Temporal_int_all_the_equations()` with the following sequence:

### 3.1 Interface & Geometry Reconstruction

**3.1.1 PLIC Reconstruction** (`VOF_reconstruct_interface`)
- Computes interface normal using the **Mixed Youngs-Centered (MYC)** method
- For each interfacial cell ($0 < F < 1$):
  1. Compute normal: $\mathbf{n}_f = \text{mycs3D}(F)$
  2. Compute plane intercept $\alpha$ such that the plane $\mathbf{n} \cdot \mathbf{x} = \alpha$ cuts the unit cube with the correct volume fraction

**3.1.2 IBM Geometry** (when `VOF_IBM` enabled)
- `Vfc_smoothing()`: Convolves solid volume fraction with a 4th-order kernel:
  $$D(x) = \frac{15}{16h}\left[\left(\frac{x}{h}\right)^4 - 2\left(\frac{x}{h}\right)^2 + 1\right], \quad |x| \leq h$$
  $$vfc_{smooth} = \frac{\sum w \cdot vfc}{\sum w}$$

- `VOF_normals_IBM()`: Computes solid surface normal from $vfc$:
  $$\mathbf{n}_s = \frac{\nabla vfc}{|\nabla vfc|}$$
  
- `VOF_normals_IBM_smooth()`: Same but using smoothed $vfc_{smooth}$

---

### 3.2 VOF Advection & Extension

**3.2.1 Geometric Flux Computation** (`VOF_set_advection`)

Uses dimension-splitting with upwind donor cell approach (Basilisk-style):
- For each face (e.g., x-direction):
  1. Determine upwind cell based on face velocity sign
  2. If upwind cell is interfacial, compute flux geometrically:
     $$c_f = \text{rectangle\_fraction}(\mathbf{n}_{plane}, \alpha, \text{lower}, \text{upper})$$
  3. Face flux: $\text{flux}_x = c_f \cdot u_{face}$

- Convective term (divergence of fluxes):
  $$\text{conv} = \frac{\partial (Fu)}{\partial x} + \frac{\partial (Fv)}{\partial y} + \frac{\partial (Fw)}{\partial z}$$

**3.2.2 RK3 Update** (`VOF_update_F`)

$$\text{RHS}^{(k)} = \gamma^{(k)} (-\text{conv}^{(k)}) + \zeta^{(k)} (-\text{conv}^{(k-1)})$$
$$F^{(k+1)} = F^{(k)} + \Delta t \cdot \text{RHS}^{(k)}$$

Values are clamped: $F \in [0, 1]$.

**3.2.3 Contact Line Extension** (`VOF_geometric_extend`)

Implements the **Liu & Ding (2015)** characteristic extension method to enforce the prescribed contact angle $\theta$ in ghost cells inside the solid region.

**Physical Setup:**
- Contact angle $\theta$ is measured through the liquid from the solid surface
- For superhydrophobic surfaces: $\theta > 90°$ (e.g., $\theta = 154°$)
- For hydrophilic surfaces: $\theta \leq 90°$

**Step 1: Identify Ghost Contact-Line Region (Eq. 18)**

A cell at $(i,j,k)$ is in the ghost region if:
1. It is inside or at the solid boundary: $vfc_{i,j,k} > 0.005$
2. It has at least one fluid neighbor ($vfc < 0.5$) — ensures we can compute gradients

**Step 2: Compute Solid Normal $\mathbf{n}_s$ (Pointing INTO Fluid)**

The IBM normal field `nx_IBM, ny_IBM, nz_IBM` is computed as $\nabla vfc / |\nabla vfc|$, which points **into the solid** (from $vfc=0$ to $vfc=1$). For the extension algorithm, we need the outward normal:
$$\mathbf{n}_s = -\frac{\nabla vfc}{|\nabla vfc|}$$

This $\mathbf{n}_s$ points from solid into fluid.

**Step 3: Compute Wall Tangent $\mathbf{t}_s$**

The tangent direction is obtained by projecting $\nabla F$ onto the wall plane:
$$\mathbf{t}_s = \frac{\nabla F - (\nabla F \cdot \mathbf{n}_s)\mathbf{n}_s}{|\nabla F - (\nabla F \cdot \mathbf{n}_s)\mathbf{n}_s|}$$

The gradient $\nabla F$ is evaluated at a **fluid neighbor** cell to ensure accurate interface orientation. This is critical because $\nabla F$ computed inside the solid may be unreliable.

**Step 4: Construct Characteristic Stencil Directions (Eq. 17)**

Two symmetric stencil directions $\mathbf{M}_1$ and $\mathbf{M}_2$ are constructed:
$$\mathbf{M}_1 = \mathbf{n}_s \sin\theta + \mathbf{t}_s \cos\theta$$
$$\mathbf{M}_2 = \mathbf{n}_s \sin\theta - \mathbf{t}_s \cos\theta$$

**Geometric interpretation:**
- The angle between $\mathbf{M}_{1,2}$ and the wall tangent is $(\theta - 90°)$
- For $\theta = 90°$: $\mathbf{M}_{1,2} = \mathbf{n}_s$ (stencils point normal to wall)
- For $\theta = 154°$: $\sin(154°) \approx 0.44$, $\cos(154°) \approx -0.90$
  - $\mathbf{M} \approx 0.44\mathbf{n}_s - 0.90\mathbf{t}_s$ (stencils tilt along the surface, into fluid)

Both directions point into the fluid region where donor values are sampled.

**Step 5: Ray Tracing to Find Interpolation Points**

For each stencil direction $\mathbf{M}$, starting from ghost cell center $\mathbf{P}$:
1. March along $\mathbf{M}$ until crossing into fluid region ($vfc < 0.5$)
2. Find intersection point $\mathbf{D}$ with a grid line connecting cell centers
3. Identify three consecutive fluid cells $A$, $B$, $C$ along the grid line for interpolation

The intersection must satisfy:
- $t > 0$ (positive distance along ray)
- Both cells adjacent to the intersected face are fluid
- Hit point lies within the cell face bounds

**Step 6: Quadratic Interpolation (Eq. 20)**

Given cells $A$, $B$, $C$ at uniform spacing $h$, and intersection point $D$ at distance $l_{AD}$ from $A$:
$$C_D = \left(\frac{l_{AD}}{h}\right)^2 \left(\frac{1}{2}C_A - C_B + \frac{1}{2}C_C\right) + \frac{l_{AD}}{h}\left(-\frac{1}{2}C_A + 2C_B - \frac{1}{2}C_C\right) + C_A$$

This is equivalent to Lagrange interpolation through points $A$, $B$, $C$.

**Step 7: Apply Extension Rule (Eq. 19)**

The extended value $C_P$ at ghost cell $P$ is determined by:
$$C_P = \begin{cases} 
\max(C_{D1}, C_{D2}) & \text{if } \theta \leq 90° \text{ (hydrophilic)} \\
\min(C_{D1}, C_{D2}) & \text{if } \theta > 90° \text{ (hydrophobic)}
\end{cases}$$

**Physical rationale:**
- Hydrophilic ($\theta \leq 90°$): Liquid spreads → maximize liquid volume fraction
- Hydrophobic ($\theta > 90°$): Liquid retracts → minimize liquid volume fraction (maximize gas)

**Step 8: Apply Constraint**

Ensure physical consistency:
$$F + vfc \leq 1$$

If violated, clip: $F = 1 - vfc$.

**Algorithm Summary:**

for each cell (i,j,k) with vfc > 0.005:

    if no fluid neighbor exists: skip
    n_s = -normalize(nx_IBM, ny_IBM, nz_IBM)  // points into fluid
    gradF = central_diff(F) at fluid neighbor
    t_s = normalize(gradF - (gradF·n_s)*n_s)

    M1 = n_s*sin(θ) + t_s*cos(θ)
    M2 = n_s*sin(θ) - t_s*cos(θ)

    C_D1 = trace_and_interpolate(origin, M1)
    C_D2 = trace_and_interpolate(origin, M2)

    if θ <= 90°:
        C_P = max(C_D1, C_D2)
    else:
        C_P = min(C_D1, C_D2)

    F[i,j,k] = clamp(C_P, 0, 1 - vfc[i,j,k])
---

**Reference:** Liu, H. & Ding, Y. (2015). "A diffuse-interface immersed-boundary method for two-dimensional simulation of flows with moving contact lines on curved substrates." *Journal of Computational Physics*, 294, 484-502.

### 3.3 Surface Tension & Contact Line Force

**3.3.1 VOF Smoothing** (`VoF_smoothing`)

Convolves $F$ with the 4th-order kernel to obtain $\tilde{F}$:
$$\tilde{F}_{i,j,k} = \frac{\sum_{stencil} w(d) \cdot F}{\sum_{stencil} w(d)}$$

This ensures $\tilde{F} \in [0,1]$ and provides smooth gradients for curvature computation.

**3.3.2 Curvature Calculation** (`curvature_patel`)

1. Compute smoothed interface normal:
   $$\hat{\mathbf{n}} = \frac{\nabla \tilde{F}}{|\nabla \tilde{F}|}$$

2. Compute curvature as negative divergence of normal:
   $$\kappa = -\nabla \cdot \hat{\mathbf{n}} = -\left(\frac{\partial n_x}{\partial x} + \frac{\partial n_y}{\partial y} + \frac{\partial n_z}{\partial z}\right)$$

3. Curvature computed only in interfacial cells ($0.005 < \tilde{F} < 0.995$); otherwise marked as `nodata`.

**3.3.3 CSF Surface Tension Force** (`Velocity_add_surfacetension_2_RHS_patel`)

Adds the **Continuum Surface Force** (Brackbill et al., 1992) to momentum RHS:
$$\mathbf{f}_{st} = \frac{2}{We} \kappa \nabla F$$

At each face (e.g., x-face at $i+1/2$):
$$\text{RHS}_u += \frac{2}{We} \cdot \kappa_{face} \cdot \frac{F_{i} - F_{i-1}}{\Delta x}$$

where:
- $We = \frac{\rho_1 U^2 L}{\sigma}$ is the Weber number
- $\kappa_{face}$ is interpolated from cell-centered values

**3.3.4 Contact Line Force (CCF)** (when `VOF_IBM` enabled)

Based on **Konstantidinis et al. (Phys. Fluids, 2024)**:

- `VOF_compute_CCF_tangents()`: Computes tangent vectors at contact line:
  $$\mathbf{t}_{int} = \frac{\mathbf{n}_s - (\mathbf{n}_f \cdot \mathbf{n}_s)\mathbf{n}_f}{|\mathbf{n}_s - (\mathbf{n}_f \cdot \mathbf{n}_s)\mathbf{n}_f|}$$
  $$\mathbf{t}_{cl} = \frac{\mathbf{n}_f - (\mathbf{n}_s \cdot \mathbf{n}_f)\mathbf{n}_s}{|\mathbf{n}_f - (\mathbf{n}_s \cdot \mathbf{n}_f)\mathbf{n}_s|}$$

- `VOF_compute_CCF_force_density()`: Computes volumetric CCF force:
  $$\mathbf{f}_{CCF} = \frac{1}{We} \cdot (\nabla F \cdot \mathbf{t}_{cl}) \cdot (\nabla vfc \cdot \mathbf{n}_s) \cdot \mathbf{t}_{int}$$

  This force is active only in contact-line cells where both $F$ and $vfc$ are interfacial.

- `Integrate_CCF_to_particle_Eulerian()`: Interpolates CCF force density to particle markers using the delta function kernel.

---

### 3.4 Gravity Force (`Velocity_add_gravity_2_RHS`)

Adds buoyancy force for two-phase flows:
$$\mathbf{f}_g = 2 \tilde{\rho} \cdot Ri \cdot \hat{\mathbf{g}}$$

where:
- $Ri = \frac{g L}{U^2}$ is the Richardson number
- $\tilde{\rho}$ is the dimensionless density field
- $\hat{\mathbf{g}}$ is the unit gravity vector

---

### 3.5 Density and Viscosity Update (`VOF_update_density_viscosity`)

Updates material properties using linear mixing laws:

**Fluid mixture** (liquid + gas):
$$\tilde{\rho}_{fluid} = F + (1-F)\frac{\rho_2}{\rho_1}$$
$$\tilde{\mu}_{fluid} = F + (1-F)\frac{\mu_2}{\mu_1}$$

**With IBM** (following O'Brien & Bussmann, JCP 2020):
Solid phase uses liquid properties for computational convenience:
$$\tilde{\rho}_{NS} = (1 - \Phi_s) \tilde{\rho}_{fluid} + \Phi_s \cdot 1$$
$$\tilde{\mu}_{NS} = (1 - \Phi_s) \tilde{\mu}_{fluid} + \Phi_s \cdot 1$$

where $\Phi_s = vfc_{smooth}$.

---

### 3.6 IBM Forcing

Applied in predictor-corrector mode within the velocity solve:

1. **Explicit velocity solve** (convective + explicit viscous terms)
2. **IBM predictor**: `Lagrangian_force(-1, ...)` adds forcing to RHS
3. **Implicit velocity solve** (implicit viscous terms)
4. **IBM corrector**: `Lagrangian_force(N_forcing_loops, ...)` iteratively corrects velocity

---

### 3.7 Pressure Projection

1. **Pressure RHS**: `Pressure_set_RHS()` computes divergence of intermediate velocity
2. **CG Solver**: `Pressure_solve_cg()` solves the variable-coefficient Poisson equation
3. **Velocity correction**: `Pressure_project_velocity()` projects onto divergence-free field

---

## 4. Summary: Chronological Method Sequence

The following table shows the **chronological order** of methods called in each RK3 substep:

| Step | Method | Function | Description |
|------|--------|----------|-------------|
| 1 | Interface Reconstruction | `VOF_reconstruct_interface` | PLIC normals & plane intercept |
| 2 | IBM Geometry (if `VOF_IBM`) | `Vfc_smoothing`, `VOF_normals_IBM` | Solid vfc smoothing & normals |
| 3 | VOF Fluxes | `VOF_set_advection` | Geometric flux computation |
| 4 | VOF Update | `VOF_update_F` | RK3 volume fraction update |
| 5 | Contact Line Extension (if `VOF_IBM`) | `VOF_geometric_extend` | Liu & Ding characteristic method |
| 6 | VOF Smoothing | `VoF_smoothing` | Smooth $F$ for curvature |
| 7 | Curvature | `curvature_patel` | $\kappa = -\nabla \cdot \hat{\mathbf{n}}$ |
| 8 | Surface Tension | `Velocity_add_surfacetension_2_RHS_patel` | CSF force to momentum RHS |
| 9 | CCF Force (if `VOF_IBM`) | `VOF_compute_CCF_tangents`, `VOF_compute_CCF_force_density` | Contact line capillary force |
| 10 | Gravity | `Velocity_add_gravity_2_RHS` | Buoyancy force |
| 11 | Density/Viscosity | `VOF_update_density_viscosity` | Update mixture properties |
| 12 | Explicit Velocity | `Velocity_solve_explicit` | Solve explicit terms |
| 13 | IBM Predictor (if `LAG_PARTICLE_RESOLVED`) | `Lagrangian_force(-1, ...)` | Initial IBM forcing |
| 14 | Implicit Velocity | `Velocity_solve` | Solve implicit viscous terms |
| 15 | IBM Corrector (if `LAG_PARTICLE_RESOLVED`) | `Lagrangian_force(N, ...)` | Iterative IBM correction |
| 16 | Pressure Solve | `Pressure_solve_cg` | CG Poisson solver |
| 17 | Velocity Projection | `Pressure_project_velocity` | Divergence-free correction |
| 18 | Particle Forces | `Lagrangian_evaluate_fluid_forces` | Collect forces on particles |
| 19 | Particle Motion | `Lagrangian_motion` | Update particle positions |

---

## 5. Summary: Physical Methods

| Component | Method | Key Reference |
|-----------|--------|---------------|
| Interface Reconstruction | PLIC with MYC normals | Scardovelli & Zaleski (1999) |
| VOF Advection | Geometric flux, dimension-split | Basilisk VOF |
| Surface Tension | CSF (Continuum Surface Force) | Brackbill et al. (1992) |
| Curvature | Smoothed Height Function variant | Popinet (2009) |
| Contact Angle Extension | Characteristic stencil method | Liu & Ding (2015) |
| Contact Line Force | CCF volumetric force | Konstantidinis et al. (2024) |
| Density/Viscosity in solid | Liquid proxy | O'Brien & Bussmann (2020) |
| IBM Forcing | Direct forcing with Lagrangian markers | Uhlmann (2005) |
| Delta Function | Roma kernel (3-point, 2nd order) | Roma et al. (1999) |

---

## 6. File Structure

| File | Purpose |
|------|---------|
| **VOF Core** | |
| `VOF_InterfaceReconstruction.c` | PLIC normals, $\alpha$ computation, `VoF_create/destroy` |
| `VOF_Advection.c` | Geometric fluxes, RK3 update, density/viscosity mixing |
| `VOF_SurfaceTension.c` | VOF smoothing, curvature, CSF force, gravity |
| **VOF-IBM Coupling** | |
| `VOF_IBM_wetting.c` | Solid normals, geometric extension, CCF tangents & force |
| **IBM / Lagrangian** | |
| `Lagrangian/Lagrangian.c` | Particle management, IBM forcing (`Lagrangian_force`, `Lagrangian_force_individual`), force collection |
| `Lagrangian/Interpolate.c` | Delta function kernel, Eul↔Lag interpolation, momentum integration, volume fraction computation |
| `Lagrangian/Particle.c` | Particle initialization, marker generation, particle lists |
| **Time Integration** | |
| `Temporal_int.c` | RK3 time loop, orchestrates all method calls |
| **Pressure** | |
| `Pressure.c` | Poisson solver, velocity projection |
| **Grid & Memory** | |
| `Cart3d.c` | Initialization, data structures, main setup |
| `Grid.c` | Grid generation, staggered grid handling |
| `Memory.c` | Memory allocation utilities |
