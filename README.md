# Diffuse-Interface Wetting FSI Extension for PARTIES

**Branch:** `VoF-MaximeJalabert`  
**Maintainer:** Maxime Jalabert (PhD Student @UCSB)  
**Status:** Active research / JCP manuscript in preparation

This branch extends the original **PARTIES** framework from a grain-resolving single-phase immersed-boundary solver to a **two-phase wetting fluid-structure interaction solver** for incompressible liquid-gas flows with fully resolved rigid bodies.

The current formulation uses a conservative **diffuse-interface (Cahn-Hilliard) phase field** with variable density and viscosity, a balanced potential-form surface-tension force, a direct-forcing **Immersed Boundary Method (IBM)** for moving solids, characteristic contact-angle enforcement on curved immersed boundaries, and a continuum capillary force whose net force and torque are fed back into the rigid-body equations.

The implementation is research-grade and under active development. The wetting machinery currently targets planar two-dimensional and axisymmetric configurations.

---

## Scope and Goals

This branch focuses on:

- Conservative Cahn-Hilliard liquid-gas interface transport with variable density and viscosity.
- Potential-form diffuse-interface surface-tension forcing, avoiding explicit VoF curvature reconstruction.
- Direct-forcing IBM coupling for fully resolved rigid bodies on a staggered Cartesian grid.
- Prescribed static contact-angle enforcement on curved immersed solids through characteristic phase-field extension into the diffuse solid.
- Continuum capillary force and torque evaluation for rigid-body feedback, alongside hydrodynamic and gravity loads.
- Validation from canonical two-phase flow benchmarks through static and dynamic wetting fluid-structure interaction.

---

## Implemented Features

### 1. Diffuse-Interface Two-Phase Solver

- Conservative transport of the liquid indicator `C_L` with a Cahn-Hilliard phase-field equation.
- Liquid, gas, and solid volume-fraction indicators `C_L`, `C_G`, and `C_S` forming the local material properties.
- Potential-form surface-tension forcing based on the liquid-gas chemical potential.
- Planar two-dimensional and axisymmetric staggered-grid discretizations advanced within a low-storage Runge-Kutta cycle.
- Validated on:
  - Axisymmetric rising-bubble convergence and mass-conservation tests.
  - Rising-bubble topology change to a toroidal bubble.
  - Single-mode Rayleigh-Taylor instability.

### 2. Immersed-Boundary Coupling for Resolved Bodies

- Direct-forcing immersed-boundary representation of fully resolved rigid particles.
- Hydrodynamic force and torque recovery from the IBM reaction force with a sharp fictitious-fluid correction.
- Rigid-body motion advanced with hydrodynamic, gravity, and capillary contributions.
- Validated on:
  - A translating/falling cylinder at `Re = 40`, including drag coefficient and wake recirculation length.

### 3. Moving Contact-Line Treatment

- Characteristic contact-angle extension of the phase field into the diffuse solid.
- Contact-angle enforcement on curved immersed boundaries in planar and axisymmetric geometries.
- Static wetting validation on:
  - Droplets on cylinders over prescribed contact angles.
  - Gravity-free droplets on spheres, including capillary and pressure-load force balance.

### 4. Continuum Capillary Force on Solids

- Eulerian continuum capillary force density evaluated from the same phase-field and diffuse-solid support.
- Net capillary force and torque passed directly to the rigid-body equations.
- Equilibrium tests recover prescribed contact angles within about `1.5%`.
- Sphere force-balance residuals and liquid-volume errors show near second-order convergence on refined grids.

### 5. Dynamic Wetting Benchmarks

- Sinking cylinder from a water surface, reproducing the contact-line migration and trajectory measured by Vella et al.
- Superhydrophobic sphere impact, reproducing the sinking-versus-bouncing bifurcation measured by Lee and Kim.
- Three-cylinder impact demonstration with multiple moving contact lines and interacting capillary cavities.

---


# PARTIES Documentation  

## Main Information
[Download and Install the code](https://github.com/metialex/PARTIES/blob/master/docs/md_files/download_and_installation.md) <br>
[Simulation setup](https://github.com/vowinckel/PARTIES/blob/master/docs/md_files/simulation_setup.md)  <br>
[Post-processing](https://github.com/metialex/PARTIES/blob/master/docs/md_files/post_proc.md) <br>

## Guidelines

[Workflow and branching](https://github.com/vowinckel/PARTIES/blob/master/docs/md_files/guidline.md) <br>
[Testcases guideline](https://github.com/metialex/PARTIES/blob/master/docs/md_files/test_case_doc.md) <br>
Writing the code - comments, description of functions <br>

## Fluid Dynamic Theory
Non-dimensionalization of code <br>
Eulerian fluid solver <br>
Lagrangian particle solver <br>
Concentration field and single-phase flow solver <br>
Immersed Boundary Method (IBM) <br>
Momo's IBM <br>
Turbulence models <br>

## Numerical 
Finite difference <br>
Time & space integration <br>
Boundary & Initial conditions <br>
Eulerian & Lagrangian solvers <br>
Paralleliztion of the code <br>

## Other
[h5ToVTK](https://github.com/metialex/h5ToVTK) - this tool converts PARTIES Eulerian and Lagrangian output .h5 files into .vtk files. <br>
[Bug reports](https://github.com/metialex/PARTIES/blob/master/docs/md_files/bug_report.md)  <br>
[Main git commands](https://github.com/metialex/PARTIES/blob/master/docs/md_files/git_commands.md) <br>
[Matlab scripts](https://github.com/metialex/PARTIES/blob/master/docs/md_files/matlab_scripts.md) <br>
[Testcases](https://github.com/metialex/PARTIES/tree/master/PARTIES/testcases/README_Testcases.md) <br>
The main description of code structure <br>
Links to other related applications - e.g. non-dimensionalization calculators, pre- and post-processing scripts. 

# Structure of the documentation
readme.md - main file with links to the main chapters

folder with different .md documents and figures <br>
folder with related papers and other (preferably published) documents/manuscripts <br>
other files (presentations, schematis ...) 
