PARTIES
Finite-differences flow solver for particulate flow using IBM.

# Background

If your are unfamiliar with Linux, please checkout this tutorial https://ryanstutorials.net/linuxtutorial/

# Howto install PARTIES

## New generic way

I added an installation script to PARTIES that manages everything. Please clone the project again, with "git clone https://github.com/tk-koellner/PARTIES.git". Then go to  PARTIES/ and type "sh  install_parties.sh".

If no errors occur, PARTIES should have all the libs and can be build by typing "make". The libs are hosted on a google drive and I download them with a Linux tool named "curl" that has to be available. If this step is broken, please download the archives manually and place them at ${BASE_DIR} https://drive.google.com/open?id=13Et77fT1G2yk9179X6uHCHANqh6D5T2M. Note that the shell variable MY_MACHINE has to be unset for this way!

Add the following line to $HOME/.bashrc

	 export PATH=$PATH:$HOME/Software/PARTIES_Libs/bin

in order to find the mpirun command.

## OLD generic way
1. With the preferred C compiler the following libraries have to be installed: **OpenMPI(or other MPI lib)**, **HDF5**, **FFTW**. The most recent lib can be found by searching the web.

2. A "make.def. $MY_MACHINE" has to to be configured with the paths to the formerly installed libs and the "mpicc" compiler. The enviroment variable $MY_MACHINE has to be present, i.e., include that to your login file ("$HOME/.bashrc")  
3. Type make to build an executable



A basic tutorial for linux is found here: "https://ryanstutorials.net/linuxtutorial/"  and here "https://www.youtube.com/watch?v=cyINirlJZzk"

A tutorial for building software: https://www.linuxvoice.com/linux-101-how-to-compile-software/

A description of how to use git can be found here https://guides.github.com/.

A tutorial about Linux shell scripts: https://bash.cyberciti.biz/guide/Main_Page.

## Stampede 2

### Environment
Here the necessary libraries are already compiled for you. They just have to be loaded with the "module" utility.

Add the following lines to your login script  .bashrc:

<ul>
	<li> export MY_MACHINE=stampede2</li>
	<li> module load fftw3 </li>
<li>module load phdf5</li>
<li>module load swr</li>
<li>module load visit</li>
<li>module load matlab</li>
 </ul>

 Source the login script or open a new terminal.

 All the information about Stampede 2 are here https://portal.tacc.utexas.edu/user-guides/stampede2.


 ### Running jobs doing post processing
 An example script to submit jobs to the queue is in the "/input_files/run_stampede2.sh". It is understood that you run  "sbatch run_stampede2.sh" in the working directory that contains "parties".


For interactive work, you could use the service here https://vis.tacc.utexas.edu/ to start a vnc session. But since the waiting time on the standard queues might be very long, it makes only sense to start a session on the "development", "skx-dev" queues here.

For longer session submit a DCV job using the script here "/share/doc/slurm/". When started, it will create a file in your home "dcvserver.out", here you find a web-adress that you can use.

*When doing large simulation, i.e. using more than 4 nodes:*

Please increase the stripe count (= number of Object Storage Targets) with

	>  lfs setstripe -c 30 $PWD

for each directory large data will be written in. New directories inherit his this propery from their parent. The status can be checked with

		> lfs getstripe .

For a description of the lustre file system see  (http://lustre.ornl.gov/lustre101-courses/content/C1/L1/LustreIntro.mp4)


A rule thumb  is using approx. 10 million gridpoints per SKX node. For simulation with a small number of timesteps more gridcells (100M) per node might be beneficial. Since we always pay for the whole node, the number of MPI task is 48 times the number of nodes. See example submission script in /inputfiles/run_stampede2.sh.

# Using Git

	1. Install git
	2. Create user account at github; ask tk-koellner to add to the online repository
	3. > git config –global  user.name $YOURNAME
	4. > git clone  https://github.com/tk-koellner/PARTIES.git
	(After your were added, you can download the project, this requires your github credential)
	5. > git  branch –a
	(see what branches are there, your are on the starred branch)
	6. > git fetch origin (download the online repo but does not change your local repo, you can look at the online repo with different tools; I use eclipse )
	7. > git checkout --track origin/$OTHER_REMOTE_BRANCH  # when the desired code is in anothe branchr
	8. > git checkout -b $OWNBRANCH  # this creates your branch you are working on
	9. > git push origin $OWNBRANCH  # create a branch on the remote
	10. > git status, git add; git commit, git push origin ...




# Configuration of the physical problem for compile time

## "./scr/Include/Boundary.h"


 Here several option to adapt the code are given.
		- Boundary condition for the flow solver
		- Boundary condition for the scalar fields
		- Switches for: Particles, Dry Particles, Concentrations fields, buoyancy terms, collision model

They are all controlled with C preprocessor commands, e.g. #define, #undefine


## In "./scr/IO/default.inp"
All physical paramters needed to run the code a pre-defined here. The majority of these parameters could be changes with the  Input script "parties.inp" during runtime; i.e., the values in default.inp might be overwritten at runtime.


## Configuration of the physical at runtime
	1. stop.inp
	2. parties.inp
	3. p_mobile.inp
	4. p_fixed.inp


###parties.inp

	1. for the scalar fields (conc_int_type)

	..1.      Lock Exchange
	..2.
	..3.
	..4.
	..5. c = sin(2*pi*x/L_x) ;
	..6. c = (1-y/Ly)+sin(2*pi/Ly*x)*1e-4

### p_mobile.inp / p_fixed.inp
In this text files all particles thta are used in the simulation are listed. In the first row the number of particles is specified. Then each row represent one particles by:

	X Y Z Radius	?? ??

#### Scripts for the initial conditions of particles


#### Random placement of particles
In "/particle/random_distribution_particle/", a matlab script is available that produces an initial distribution of particles. Input parameters are the median radius, the number of particles the space they occupy, etc.

#### Hexagonal packing


## Running parties

### Local machines


### Knot cluster



For a description of the lustre file system see  (http://lustre.ornl.gov/lustre101-courses/content/C1/L1/LustreIntro.mp4)

# Visualization

# Visualization

## XDMF Writer
Produces a file for Vist and Paraview.

## Matlab scripts




## What are output files


#Prgramming

## Add measures for analyzing

To include new  outputs do the following:

 1. create new switch in parameters to control if  calculation and output is performed
 2. add this parameter to default.inp  
 3. allocate memory in Statistics2d_create() for the new variable also destroy memory in Statistics2d_destroy()
 4. add calculation to Statistics2d_computeStatistics()
 5.  do the  output in void Statistics2d_writeStatistics2dToH5()  
 6. add the so defined variable to the manual
 7. check if calculations are correct!!





## TO DO

- output energy and momentum balance for particles (documentation)
- compute energy budget with buoyancy (available pot energy)
