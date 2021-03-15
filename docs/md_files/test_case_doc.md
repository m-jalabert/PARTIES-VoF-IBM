# Testcases guideline
Testcases are required to keep the consistency of the code and avoid undesireable changes in the base parts of the code 
which are responsible for output results (such as numerical solvers, discretization schemes and others).

The testcases are used in case of a pull request which includes changed or added numerical functionality.
Before merging your branch with a master branch you have to perform all the testcases simulations to check 
if your additions to the code affected the results of the testcases simulations.

It is also possible that the results could be different than before. In this case you have to clearly state 
the reason and provide a strong motivation for your changes.

In case if your updated code provides an additional functionality, it is highly encouraged to create a new test case,
which would represent a simulation where your new functionality is tested. Below you can find a section which describes the exact 
procedure how a new test case should be added (section __Add a testcase__).


## Testcases structure
The directory `<$PARTIES_home>/testcases` is assigned for the testcases related files and folders.

The following folders and files are located in `<$PARTIES_home>/testcases`:
 - The folders which correspond to a certain testcase (e.g. Poiseuille_Testcase(link)).
 - The main script file which runs and checks all testcases runtestcases.sh(link).
 - The folder which contains validated results of testcases Testcase_Results(link).

The folder with a certain test case contains the following files and folders:
 - The setup files which are used to perform a certain testcase (e.g. parties.inp(link), Boundary.h(link), p_fixed.inp(link)).
 - The individual script file which is executed to run and check a corresponding testcase (e.g. Poiseuille.sh(link)).
 - Files which contain the data for initial validation of a testcase. These are usually the files related to the output data from anlitical solution or experimental results.


## Run testcases
### All testcases
To run all test cases, go to the folder testcase(add link), open terminal and run the script which includes all test cases:
```
$ cd <$PARTIES_home>/testcases/
$ ./runtestcases.sh > logfile.log
```
The script `runtestcases.sh` includes the list of all test cases and executes them in a sequential manner.
The status of each test case is written in the logfile `logfile.log`.
Therefore if the status of every test case is "passed" the code could be submitted for a pull request.
In case if some of test cases did not passed you might think about the reasons of it.
If you have a certain explanation for it, for example you implemented a new method for solving PDE, try to compare your results with analitical or experimental results.
In case if you dont have an explanation why results are different - try to double check your code and find out why the results are different.

__Note :__ the intermediate results for each test case are deleted after the passing criteria is checked. 

### Certain testcase 
If you want to analize the simulation data for a certain testcase in more details, run a script for this testcase (e.g. Poiseuille.sh(link)).
Each individual script uses the simulation setup files from corresponding folder and compiles the PARTIES code based on them.
Simulations starts after the compilation the code. Finally, the ouput data is compared with the results file from the same simulation.
The simulation results are stored in the work folder which is automatically created in `<$PARTIES_home>/testcases/ (check it)`.


## Passing criteria
(in progress)

what data is compared (velocity profiles for example)
Simulation data compared with previous simulation data
Simulation is also compared with analitical solution/experiments but it is not a passing check criteria (but maybe initial criteria to add a test case)

## Add a testcase
(in progress)

# Ralphs readme.md
## PARTIES Testcase for standard flows
This shell script runs the PARTIES CFD solver for a standard flow and compares the obtained numerical solution against a verified numerical solution 
of the same flow type. 
The initial verification is done against the analytical solution of the fluid flow.
Eg: Poiseuille flow, Couette flow


### Prerequisites
Valid PARTIES installation.

Files:
<flowname>_testcase.sh
Boundary_<flowname>.h
p_fixed.inp
p_mobile.inp
parties.inp
stop.inp
xdmfWriter.inp
l2norm.c
numerical_velo_verified_<flowname>.dat


### Initial Setup
Analytical solution:
<flowname>_analytical_soln.m to be executed beforehand, to obtain the 'analytical_velo_<flowname>.dat' file
This is used to verify the numerical soultion, obtained through the PARTIES solver, in order to acquired the numerical_velo_verified_<flowname>.dat file.


### Usage
The files mentioned above are present in <$PARTIES_home>/testcases/<flowname>_Testcase directory.

In order to run multiple testcases in series, the bash "script runtestcases.sh" is present.
Located in <$PARTIES_home>/testcases/ directory, it executes the testcases present in a sequential manner. 

provide execution rights to the bash script, if needed.
$ chmod 777 <flowname>_testcase.sh

Result is the L2 Norm between the verified numerical velocities and numerical velocities along the y-direction of the control volume.


### Output
For passed condition, when the L2 norm is with the allowed tolerance, the numerical_velo_verified_<flowname>.dat is updated, so as to be the one from the lastest simulation.
For an unfavourable test result, the numerical_velo_verified_<flowname>.dat is not updated.
The ouput is place in the <$PARTIES_home>/testcases/Testcase_Results directory.


### Contribution
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests and README.md as appropriate.


### License
Add type of licens
