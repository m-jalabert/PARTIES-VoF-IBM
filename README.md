# Volume of Fluid (VoF) Multiphase Extension for PARTIES

**Branch:** `VoF-MaximeJalabert`  
**Maintainer:** Maxime Jalabert (PhD Student @UCSB)  
**Status:** Active research / beta

This branch extends the original **PARTIES** framework from a single-phase particle-laden flow solver to a **multiphase Volume of Fluid - Piecewise Linear Interface Construction (VoF–PLIC) + Immersed Boundary Method (IBM)** solver capable of handling free surfaces, moving contact lines, and fluid–particles interaction.

The implementation is research-grade (under active development), but already supports a set of canonical multiphase benchmarks and contact-line problems used in my PhD work.

---

## Scope and Goals

This branch focuses on:

- Geometric **VoF–PLIC** interface tracking for multiphase flows.
- **Surface tension** and curvature modeling via a Continuum Surface Force (CSF) approach.
- Coupling VoF with the existing **Immersed Boundary Method** for moving solids.
- Validated test cases (Standard bubble advection, static droplet (Laplace pressure test), rising bubble, equilibrium droplet on a spherical surface) as a basis for more complex physics (moving contact line).

---


## Implemented Features

### 1. Interface Tracking (VoF–PLIC)

- Piecewise Linear Interface Construction (PLIC) for the volume fraction field `F`.
- Directional sweeps with conservative flux computation to maintain global mass conservation.
- Interface sharpening consistent with the underlying discretization.
- Validated on:
  - **Standard Bubble Advection Test**

### 2. Surface Tension and Curvature

- **CSF formulation** (Continuum Surface Force) using reconstructed interface normals.
- Curvature computed from VoF-derived normals and fed back as a body force in the momentum equation.
- Validated on:
  - **Static droplet**:
    - Laplace pressure jump vs. theoretical prediction.
    - Spurious current levels used as a quality metric.
  - **Rising bubble**:
    - Terminal rise velocity and bubble shape compared against reference data.

### 3. Immersed Boundary Method (IBM) Coupling for static bodies

- Coupling of VoF with the existing IBM framework for static bodies.
- Volume of Fluid extension via an advection PDE and contact angle enforcement at the contact line via rotation of the smoothed normals.
- Validated on:
    - **Equilibrium droplet on a spherical surface**, for contact angles from 30° to 150°.

### 4. Moving Contact Line — In progress

- Implementation of the Continuum Capillary Force (CCF) for the solid equations.
- Replacement of the older **PDE-based extension** by a **geometric extension** strategy for `F` in solid cells.
- Ongoing work toward:
  - Robust handling of moving contact lines on spherical immersed geometries.
- Test case to be handled:
  - **Impact of a superhydrophobic sphere**
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
