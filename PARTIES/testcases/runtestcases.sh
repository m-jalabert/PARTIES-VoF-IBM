#!/bin/bash
# Shell script to run selected test flows with the PARTIES fluid solver.
# v1
# Ralph George: r.george@tu-braunschweig.de

# Current directory will be <PARTIES_home>/testcases/
echo '----------**Starting PARTIES Testcases**----------'
test_home_path=`pwd`

# Path to executables
POISEUILLE_FLOW="$test_home_path/Poiseuille_Testcase/poiseuille_flow_testcase.sh"
COUETTE_FLOW="$test_home_path/Couette_Testcase/couette_flow_testcase.sh"

# Input no of processors
read -p 'Enter number of processors: ' nproc

###########################################################################
#                          Testcase Execution                             #
###########################################################################

# Comment out unrequired flows
printf '\n\n'
cd Poiseuille_Testcase/
. "$POISEUILLE_FLOW" $nproc

printf '\n\n'
cd Couette_Testcase/
. "$COUETTE_FLOW" $nproc

###########################################################################
#                                Clean up                                 #
###########################################################################
cd $test_home_path
rm -rf Poiseuille_run
rm -rf Couette_run
