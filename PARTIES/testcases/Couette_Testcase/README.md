###########################################################################
#                   PARTIES Testcase for standard flows                   #
###########################################################################
The shell script runs the PARTIES solver and compares it against a standard analytically solved fluid flow.
Eg: Poiseuilli flow, Couette flow


###########################################################################
#                              Prerequisites                              #
###########################################################################
Valid PARTIES installation.

Files:
<flowname>_testcase.sh
l2norm.c
p_fixed.inp
p_mobile.inp
parties.inp
stop.inp
xdmfWriter.inp


###########################################################################
#                                  Usage                                  #
###########################################################################
Folder 'Poiseuilli_Testcase', with the files mentioned above to be placed in testcases/<Flow Type> in working directory.

Analytical solution:
<flowname>_analytical_soln.m to be executed beforehand, to obtain the 'analytical_velo.dat' file

provide execution rights to the bash script, if needed.
$ chmod 447 <flowname>_testcase.sh

Result is the L2 Norm between the analytical and numerical velocities along the y-direction of the control volume.


###########################################################################
#                             Contribution                                #
###########################################################################
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
Please make sure to update tests and README.md as appropriate.


###########################################################################
#                                License                                  #
###########################################################################
# Add type of license
