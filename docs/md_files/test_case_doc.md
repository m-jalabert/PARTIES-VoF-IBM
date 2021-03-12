# Testcases guideline
Testcases are required to keep the consistency of the code and avoid undesireable changes in the base parts of the PARTIES code 
which are responsible for numerical solvers.

The testcases are used in case if you changed a code and want to merge it with the master branch.
In this case before merging your code you have to perform all the testcases simulations to check 
if your additions to the code affected the results of the testcases simulations.

It is also possible that the results could be different than before, but in this case you have to clearly state 
the reason and and provide a strong motivation for your changes.

In case if your updated code provides an additional functionality, it is highly encouraged to create a new test case,
which would represent a simulation where your new functionality is tested. Below you can find a section which describes the exact 
procedure how a new test case should be added.


## Testcases structure

The directory `<$PARTIES_home>/testcases`  contains all test cases and related files.
The script `testcases.sh` runs all test cases and provides their status.
Inside each testcase there are the setupfiles (e.g. parties.inp, Boundary.h) which are used to perform a certain testcase.
In addition, in each test case folder there is a script which runs and checks this certain test case.

## Passing criteria


## Run testcases

To run all the test cases, go to the folder testcase(add link), open terminal and run the script which includes all test cases:
```
$./testcases.sh > logfile.log
```
The script `testcases.sh` includes the list of all test cases and executes them in a sequential manner.
The intermediate results for each test case are deleted after the passing criteria is checked.

The status of each test case is written in the logfile `logfile.log`.
Therefore if the status of every test case is "passed" the code could be submitted for a pull request.
In case if some of test cases did not passed you might think about the reasons of it.
Please note, when you run 

## Add a testcase


# PARTIES Testcase for standard flows
This shell script runs the PARTIES CFD solver for a standard flow and compares the obtained numerical solution against a verified numerical solution 
of the same flow type. 
The initial verification is done against the analytical solution of the fluid flow.
Eg: Poiseuille flow, Couette flow


## Prerequisites
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


## Initial Setup
Analytical solution:
<flowname>_analytical_soln.m to be executed beforehand, to obtain the 'analytical_velo_<flowname>.dat' file
This is used to verify the numerical soultion, obtained through the PARTIES solver, in order to acquired the numerical_velo_verified_<flowname>.dat file.


## Usage
The files mentioned above are present in <$PARTIES_home>/testcases/<flowname>_Testcase directory.

In order to run multiple testcases in series, the bash "script runtestcases.sh" is present.
Located in <$PARTIES_home>/testcases/ directory, it executes the testcases present in a sequential manner. 

provide execution rights to the bash script, if needed.
$ chmod 777 <flowname>_testcase.sh

Result is the L2 Norm between the verified numerical velocities and numerical velocities along the y-direction of the control volume.


## Output
For passed condition, when the L2 norm is with the allowed tolerance, the numerical_velo_verified_<flowname>.dat is updated, so as to be the one from the lastest simulation.
For an unfavourable test result, the numerical_velo_verified_<flowname>.dat is not updated.
The ouput is place in the <$PARTIES_home>/testcases/Testcase_Results directory.


### Contribution
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests and README.md as appropriate.


### License
Add type of licens
