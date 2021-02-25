# Simulation setup

The setup for the simulation is done using the input files. Each input file feeds
values to the code, dealing with specific input parameters. To be noted: For a
single run, all input files are to be present. Best practice: Once the code has
been compiled, the executable, along with the input files are copied to an active
folder(run), where the output(results data files) are generated.
The default values for each of the input parameters are present in the
default.inp, present at [src/Eulerian/IO/default.inp][1]

## [parties.inp][2]

- **Geometry**: The dimensions of the control volume is setup, as per required,
using Cartesian coordinates.

- **Grid**: The values for the grid is populated to keep the number of elemental
nodes per unit volume constant and a multiple of 10, along each of the 3
coordinate axes.

- **Flow**: The Reynold’s number and the target ubulk velocity is populated.

- **Simulation**: Parameters related to time are setup here. See subsection **Simulation parameter** for infor-
mation on the various sub parameters.

- **Particle**: ADD DATA HERE
- **Concentration**: ADD DATA HERE
- **Lock**: ADD DATA HERE
- **Output**: ADD DATA HERE

### Simulation parameter

- **time_max** is max simulation time in seconds.

- **output_time_interval** On screen time interval of output.

- **cfl** the Courant–Friedrichs–Lewy condition number.

- **max_dt** max allowed time step with cfl active.

- **default_dt** default time step with cfl inactive and starting time step with
cfl active.

- **constant_dt** boolean switch, cfl active or inactive.

- **resume** defines start time of simulation

## [stop.inp][3]

Input file, which is used as an emergency switch, to stop a run, as and if required.
This file is check before the start of every main loop iteration of the flow solver.

0 --> continue <br>
1 --> stop

## [p_fixed.inp][4] and [p_mobile.inp][5]
These input files describe the initial conditions of fixed and mobile particles.

Both files are structured in the following way:
```
n
x1 y1 z1 R1
x2 y2 z2 R2
...
xn yn zn Rn
```
where `n` is a total number of particle, `x y z` correspond to the cartesian coordinates of a particle and `R` represents the radius of a particle <br>

The total number of particles `n` defines how many lines must be red in the `p_fixed.inp` and `p_mobile.inp` files. This means that if `n` is more than the number of lines which contain the particle information, then code will throw an error. If `n` is less than the number of lines which contain the particle information, only `n` particles will be considered.


## [xdmfWriter.inp][6]
ADD DATA HERE


## [Boundary.h][7]

The Boundary.h file is the main setup file for the control volume. It contains
predefined headers for selection of required boundary conditions, selection of
flow solver, implementing particles into the fluid and their concentration and
the selection turbulence models for flows. These options are to be selected as
per required for a given computational run.

### Boundary Conditions

The setup for the control volume is primarily done with the selection of the type
of boundary condition which is required. This is further broken down to the
specific conditions required along the 3 coordinate axes as well as the four/six
faces of the control volume in question.


### Implemented Boundary Conditions

These sets of boundary conditions are predefined in the code and are for the
various faces of the control volume. An option to manually setup the individual
conditions for faces of the control volume is provided as well.
1. PERIODIC NOSLIP BOX
2. PERIODIC FREESLIP BOX
3. XPERIODIC FREESLIP BOX
4. PERIODIC SHEAR FLOW
5. PERIODIC DOUBLE SHEAR FLOW
6. PERIODIC DOUBLE SHEAR FLOW
7. NOSLIP BOX
8. FREESLIP DUCT
9. NOSLIP DUCT
10. GRAVITY CURRENT
11. GRAVITY CURRENT PERIODIC
12. TRIPLE PERIODIC
13. TRIPLE PERIODIC SHEAR
14. LEFT RIGHT INFLOW TOP OUTFLOW

Reminder:
- No-slip boundary condition: At the interface between a moving fluid
and a stationary wall, both the normal and tangential components of the
fluid velocity field are equal to zero.

- Free-slip boundary condition: At the interface between a moving fluid
and a stationary wall, the normal component of the fluid velocity field is
equal to zero, but the tangential component is unrestricted. This condition
is also know as the no-penetration condition.

- Periodic boundary conditions: The flows across two opposite planes
of the control volume model are identical



[1]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/src/IO/default.inp
[2]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/parties.inp
[3]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/stop.inp
[4]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/p_fixed.inp
[5]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/p_mobile.inp
[6]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/xdmfWriter.inp
[7]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/src/Include/Boundary.h

