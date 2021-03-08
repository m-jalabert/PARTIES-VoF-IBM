# Post-processing

The output data files from simulation are the .h5 files. They are divided into two different groups: Eulerian (fluid and concentration) and Lagrangian (particles) data.
Additionally, files with the extension .dat are used for output of table-stuructured data (e.g. mobile.dat).

__Important notice__: as far as low storage Runge-Kutta scheme is used as a time integration scheme, all output variables (excluding some variables which are created if the flag "POST_PROCESS" is defined) correspond to the variables of the third iteration of RK scheme. 

The main categories of output data are described below:

## Data.h5
This set of files contain the main information about the Eulerian field of a simulation domain.<br>
The structure of the file is represented below and can be edited [here](https://github.com/metialex/PARTIES/blob/master/docs/other/Data_h5_structure.drawio) using draw.io:

![](https://github.com/metialex/PARTIES/blob/master/docs/figures/Data_h5_structure.png)

Here,<br>
`NX NY NZ` - is the number of cells in x,y and z direction,<br>
`xc yc zc` - are the arrays with coordinates of cell centers,<br>
`xu yv zw` - are the arrays with coordinates of face centers where u,v and w components are stored respectively (see [staggered grids](https://www.cfd-online.com/Wiki/Staggered_grid)),<br>
`time` - contains a double variable corresponding to a current time,<br>
`p u v w` - are the sets of three dimensional arrays contatining the data of pressure and velocity components.


## Particle.h5
This set of files contains the main information about the Lagrangian data - the main data corresponding to each particle.<br>
The structure of the file is represented below and can be edited [here](https://github.com/metialex/PARTIES/blob/master/docs/other/Particle_h5_structure.drawio) using draw.io:

![](https://github.com/metialex/PARTIES/blob/master/docs/figures/Particle_h5_structure.png)

Here,<br>
`periodic` - is a flag variable, which indicates active periodic BC for each direction,<br>
`xmin xmax ymin ...` - are variables defining the global size of the domain,<br>
`time` - contains a double variable corresponding to a current time,<br>

Data sets `fixed` and `mobile` contains the variables for each particle.<br>
There is a big list of variables which includes `F` - force vector, `X` position vector, `U` velocity vector and other variables.

## mobile.dat
This file contains the positions and velocities of mobile particle over the time.
In contrast to `Particle.h5` output files, the _mobile.dat_ file contains the data from all time steps independently of output interval variable.

