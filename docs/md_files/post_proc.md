# Post-processing

The output data files from simulaation are the .h5 files. They are divided into two different groups: Eulerian and Lagrangian data.
Additionally, files with the extension .dat are used for output of table-stuructured data (e.g. mobile.dat)
The main categories of output data are described below:

# Data.h5
This set of files contain the main information about the Eulerian field of a simulation domain.
The structure of the file is represented here: (link)

Here
`NX NY NZ` - number of cellx in x,y and z direction
`xc yc zc` - arrays with cell center coordinates
`xu yv zw` - arrays with coordinates of face centers where u,v and w components are stored respectively (see [staggered grids](https://www.cfd-online.com/Wiki/Staggered_grid))
`time` - contains a double variable corresponding to a current time
`p u`


# Particle.h5

# mobile.dat
