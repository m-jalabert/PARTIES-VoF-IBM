# PARTIES Testcase for standard flows
The shell script <$PARTIES_home>/testcases/runtestcases.sh executes the PARTIES CFD solver for a standard flow and compares the obtained numerical solution against a verified numerical solution of the same flow type. 
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

Additional file:
MATLAB script for the analtical solution is present, along with the output.

## Initial Setup
Analytical solution:
<flowname>_analytical_soln.m to be executed beforehand, to obtain the 'analytical_velo_<flowname>.dat' file
This is used to verify the numerical soultion, obtained through the PARTIES solver, in order to acquired the numerical_velo_verified_<flowname>.dat file.

## Usage
The files mentioned above are present in <$PARTIES_home>/testcases/<flowname>_Testcase directory.

This script executes the testcases present in a sequential manner. 

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
Add type of license
