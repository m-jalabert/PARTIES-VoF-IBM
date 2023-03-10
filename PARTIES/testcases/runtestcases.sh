#!/bin/bash
# Shell script to run selected test flows with the PARTIES fluid solver.
# v1
# Ralph George: r.george@tu-braunschweig.de

#Define status print function
print_status() {
    if [[ $1 -eq 0 ]]; then tput setaf 2; printf "\t PASSED"
    else tput setaf 1; printf "\t FAILED"; fi
    tput sgr0
}

# Current directory will be <PARTIES_home>/testcases/
echo '----------**Starting PARTIES Testcases**----------'
test_home_path=`pwd`
LOGFILE="$test_home_path/main.log"

# Path to executables
POISEUILLE_FLOW="$test_home_path/Poiseuille_Testcase/poiseuille_flow_testcase.sh"
COUETTE_FLOW="$test_home_path/Couette_Testcase/couette_flow_testcase.sh"
MORDANT_TESTCASE="$test_home_path/Mordant_Testcase/mordant_testcase.sh"
TEN_CATE_TESTCASE="$test_home_path/ten_Cate_testcase/ten_Cate_testcase.sh"
COLLISION_TESTCASE_ST_3="$test_home_path/Gondret_Testcase/test_St_3.sh"
COLLISION_TESTCASE_ST_120="$test_home_path/Gondret_Testcase/test_St_120.sh"
OSCILLATION_DOMAIN="$test_home_path/Oscillation_Testcases/01_Osc_Domain/Osc_Dom_testcase.sh"
OSCILLATION_PARTICLE="$test_home_path/Oscillation_Testcases/02_Osc_Particle/Osc_Par_testcase.sh"
RHE_VOL_IMPOSED="$test_home_path/Rheology_Testcases/Volume_Imposed_Testcase/volume_imposed_testcase.sh"
RHE_PRES_IMPOSED="$test_home_path/Rheology_Testcases/Pressure_Imposed_Testcase/Pressure_imposed_testcase.sh"



# Input no of processors
read -p 'Enter number of processors: ' nproc

###########################################################################
#                          Testcase Execution                             #
###########################################################################

# Comment out unrequired flows
printf '\nPoiseuille flow\t'
cd Poiseuille_Testcase/
. "$POISEUILLE_FLOW" $nproc > $LOGFILE
print_status $?

printf '\nCouette flow\t'
cd Couette_Testcase/
. "$COUETTE_FLOW" $nproc >>  $LOGFILE
print_status $?

printf '\nMordant test\t'
cd Mordant_Testcase/
. "$MORDANT_TESTCASE" $nproc >>  $LOGFILE
print_status $?

printf '\nten Cate test\t'
cd ten_Cate_testcase/
. "$TEN_CATE_TESTCASE" $nproc >> $LOGFILE
print_status $?

printf '\nCollision test St=3'
cd Gondret_Testcase/
. "$COLLISION_TESTCASE_ST_3" $nproc >> $LOGFILE
print_status $?

printf '\nCollision test St=120' 
cd Gondret_Testcase/
. "$COLLISION_TESTCASE_ST_120" $nproc >> $LOGFILE
print_status $?

printf '\nOscillation-Domain' 
cd Oscillation_Testcases/01_Osc_Domain/
. "$OSCILLATION_DOMAIN" $nproc >> $LOGFILE
print_status $?

printf '\nOscillation-Particle' 
cd Oscillation_Testcases/02_Osc_Particle/
. "$OSCILLATION_PARTICLE" $nproc >> $LOGFILE
print_status $?

printf '\nRheology- Volume Imposed\t' 
cd Rheology_Testcases/Volume_Imposed_Testcase/
. "$RHE_VOL_IMPOSED" $nproc >> $LOGFILE
print_status $?

printf '\nRheology- Pressure Imposed\t' 
cd Rheology_Testcases/Pressure_Imposed_Testcase/
. "$RHE_PRES_IMPOSED" $nproc >> $LOGFILE
print_status $?

printf '\n'
###########################################################################
#                                Clean up                                 #
###########################################################################
cd $test_home_path
rm -rf Poiseuille_run
rm -rf Couette_run
rm -rf Mordant_run
rm -rf ten_Cate_run
rm -rf test_St_3_run
rm -rf test_St_120_run
rm -rf Osc_Dom_run
rm -rf Osc_Par_run
rm -rf Volume_Imposed_run
rm -rf Pressure_Imposed_run

rm $LOGFILE
