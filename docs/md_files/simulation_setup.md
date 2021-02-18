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

## [p_fixed.inp][4]
ADD DATA HERE

## [p_mobile.inp][5]
ADD DATA HERE

## [xdmfWriter.inp][6]
ADD DATA HERE


[1]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/src/IO/default.inp
[2]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/parties.inp
[3]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/stop.inp
[4]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/p_fixed.inp
[5]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/p_mobile.inp
[6]: https://github.com/vowinckel/PARTIES/blob/master/PARTIES/xdmfWriter.inp
